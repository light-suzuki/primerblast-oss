"""Progress must preserve evidence and reusable reads must preserve bases."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import threading

import pytest

from primerblast_oss import progress, workflows
from primerblast_oss.genome import Genome, revcomp
from primerblast_oss.webapp import server


def test_stage_times_and_context_do_not_change_engine_results(monkeypatch):
    now = [10.0]
    monkeypatch.setattr(progress.time, 'monotonic', lambda: now[0])
    events = []
    with progress.observe(events.append) as recorder:
        progress.report('blast', database='synthetic', hypothesis=1)
        now[0] = 12.0
        progress.report('realign', 0, 4, primer='F')
        now[0] = 15.0
        progress.report('realign', 4, 4)
        timing = recorder.finish()
    assert timing == {'total_seconds': 5.0, 'stage_seconds': {'blast': 2.0, 'realign': 3.0}}
    assert events[-2]['progress']['completed'] == 4
    assert events[-2]['progress']['database'] == 'synthetic'
    assert not progress.active()


def test_observers_are_isolated_between_worker_threads():
    barrier = threading.Barrier(2)
    def run(name):
        events = []
        with progress.observe(events.append):
            barrier.wait()
            progress.report('blast', database=name)
        return events
    with ThreadPoolExecutor(max_workers=2) as pool:
        a = pool.submit(run, 'a')
        b = pool.submit(run, 'b')
        assert a.result(timeout=5)[0]['progress']['database'] == 'a'
        assert b.result(timeout=5)[0]['progress']['database'] == 'b'


def test_job_partial_snapshot_and_error_timing(monkeypatch):
    jobs = server.JobManager()
    recorded = []
    def handler(params):
        payload = {'mode': 'check', 'partial': True, 'results': [{'db': 'synthetic'}]}
        progress.report('blast')
        progress.partial(payload)
        payload['results'][0]['db'] = 'mutated'
        recorded.append(jobs.get(next(iter(jobs._jobs))))
        raise ValueError('second database failed')
    monkeypatch.setitem(workflows.HANDLERS, 'check', handler)
    monkeypatch.setattr(server.threading.Thread, 'start', lambda self: self.run())
    job = jobs.get(jobs.submit('check', {}))
    assert recorded[0]['partial_result']['results'][0]['db'] == 'synthetic'
    assert recorded[0]['partial_revision'] == 1
    assert job['status'] == 'error'
    assert 'blast' in job['timing']['stage_seconds']
    assert job['timing']['total_seconds'] >= 0
    assert '_progress_at' not in job


def fixture_genome(tmp_path):
    path = tmp_path / 'synthetic.fa'
    path.write_bytes(b'>chr1\r\nACGT\r\nTGCA\r\nAAAA\r\n')
    Path(str(path) + '.fai').write_text('chr1\t12\t7\t4\t6\n')
    return Genome(str(path))


def test_nested_session_reuses_handle_and_closes_on_failure(tmp_path, monkeypatch):
    genome = fixture_genome(tmp_path)
    expected = genome.fetch('chr1', 3, 10)
    import builtins
    original_open = builtins.open
    handles = []
    def counted_open(path, *args, **kwargs):
        handle = original_open(path, *args, **kwargs)
        if str(path) == genome.fasta:
            handles.append(handle)
        return handle
    monkeypatch.setattr(builtins, 'open', counted_open)
    with pytest.raises(ValueError, match='test failure'):
        with genome.read_session():
            assert genome.fetch('chr1', 3, 10) == expected
            with genome.read_session():
                assert genome.fetch('chr1', 3, 10, '-') == revcomp(expected)
            assert genome.fetch('chr1', 3, 10) == expected
            raise ValueError('test failure')
    assert len(handles) == 1 and handles[0].closed
    assert genome.fetch('chr1', 3, 10) == expected
    assert len(handles) == 2


def test_shared_genome_sessions_have_separate_thread_handles(tmp_path):
    genome = fixture_genome(tmp_path)
    barrier = threading.Barrier(2)
    def run(start, end, strand):
        expected = genome.fetch('chr1', start, end, strand)
        with genome.read_session():
            barrier.wait()
            for _ in range(100):
                assert genome.fetch('chr1', start, end, strand) == expected
    with ThreadPoolExecutor(max_workers=2) as pool:
        futures = [pool.submit(run, 1, 7, '+'), pool.submit(run, 4, 12, '-')]
        for future in futures:
            future.result(timeout=5)


def test_completed_databases_are_published_before_next_search(monkeypatch):
    def search(primers, database, **kwargs):
        return {'db': database, 'sites_per_primer': {}, 'n_products': 0,
                'products': [], 'search_complete': True, 'search_completeness': 'complete'}
    monkeypatch.setattr(workflows, 'in_silico_pcr', search)
    monkeypatch.setattr(workflows, '_associated_genomes', lambda params: {})
    params = {'db': ['synthetic_a', 'synthetic_b'], 'forward': 'ACGA',
              'reverse': 'CGTA', 'input_orientation': 'auto'}
    events = []
    with progress.observe(events.append):
        observed = workflows._run_check(params)
    snapshots = [event['partial_result'] for event in events if 'partial_result' in event]
    assert len(snapshots) == 8
    assert [snapshot['completed_units'] for snapshot in snapshots] == list(range(1, 9))
    assert all(snapshot['partial'] and snapshot['total_units'] == 8 for snapshot in snapshots)
    assert len(snapshots[0]['results']) == 1
    assert [row['db'] for row in observed['results'][:2]] == ['synthetic_a', 'synthetic_b']
    assert observed == workflows._run_check(params)
    assert 'partial' not in observed


def test_alignment_cache_is_bounded_and_keys_include_both_sequences():
    from primerblast_oss.specificity import _fitting_align_primer
    _fitting_align_primer.cache_clear()
    cases = [('ACGA', 'TTACGATT'), ('ACGA', 'TTACGCTT'), ('CGTA', 'TTACGATT')]
    for primer, target in cases:
        expected = _fitting_align_primer.__wrapped__(primer, target)
        assert _fitting_align_primer(primer, target) == expected
        assert _fitting_align_primer(primer, target) == expected
    info = _fitting_align_primer.cache_info()
    assert info.hits == 3 and info.misses == 3 and info.maxsize == 4096


def test_product_export_uses_resolved_genome_without_blastdbcmd(tmp_path, monkeypatch):
    from primerblast_oss.specificity import Amplicon
    from primerblast_oss import sequence_tools
    genome = fixture_genome(tmp_path)
    monkeypatch.setattr(workflows, '_associated_genomes', lambda params: {})
    monkeypatch.setattr('primerblast_oss.pipeline.resolve_genome_for_database',
                        lambda *args, **kwargs: (genome, 'db_associated'))
    monkeypatch.setattr(sequence_tools.shutil, 'which', lambda name: None)
    monkeypatch.setattr(workflows, 'in_silico_pcr', lambda primers, db, **kwargs: {
        'db': db, 'sites_per_primer': {'F': 1, 'R': 1}, 'n_products': 1,
        'products': [Amplicon('chr1', 1, 8, 8, 'F', 'R', 0, 0)],
        'search_complete': True, 'search_completeness': 'complete'})
    result = workflows._run_check({'db': ['synthetic'], 'forward': 'ACGT', 'reverse': 'TGCA'})
    product = result['results'][0]['products'][0]
    assert product['sequence_source'] == 'genome_fasta'
    assert product['sequence'] == 'ACGTTGCA'
    assert len(product['sequence']) == product['size']
