"""Linux real-genome batch/legacy equality, timing and sampled process-tree RSS."""
import argparse
from dataclasses import asdict
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from primerblast_oss.genome import Genome
from primerblast_oss.specificity import SpecParams, priming_sites_with_stats, screen_primers_with_stats


def primers_for(genome, count):
    chroms = [chrom for chrom in genome.chroms() if chrom in ('1', '2', '3', '4', '5')]
    if len(chroms) != 5:
        raise ValueError('This benchmark requires the five TAIR10 nuclear chromosomes')
    primers = {}
    for index in range(count):
        chrom = chroms[index % 5]
        length = genome.length(chrom)
        pos = 1000 + (index // 5 + 1) * (length - 3000) // (math.ceil(count / 5) + 1)
        sequence = genome.fetch(chrom, pos, pos + 19, '+').upper()
        while set(sequence) - set('ACGT'):
            pos += 100
            if pos + 19 > length:
                raise ValueError('Could not select an unambiguous primer')
            sequence = genome.fetch(chrom, pos, pos + 19, '+').upper()
        primers['P%05d' % index] = sequence
    return primers


def worker(args):
    genome = Genome(args.genome)
    primers = primers_for(genome, args.count)
    params = SpecParams(num_threads=args.threads)
    start = time.perf_counter()
    if args.worker == 'batch':
        sites, stats = screen_primers_with_stats(primers, args.db, params, 'blastn', genome)
    else:
        sites, stats = [], {}
        for name, sequence in primers.items():
            found, metadata = priming_sites_with_stats(sequence, name, args.db, params, 'blastn', genome)
            sites.extend(found)
            stats[name] = metadata
    elapsed = time.perf_counter() - start
    canonical = json.dumps({'sites': [asdict(site) for site in sites],
                            'stats': {name: asdict(value) for name, value in stats.items()}}, sort_keys=True)
    return {'seconds': elapsed, 'n_sites': len(sites),
            'evidence_sha256': hashlib.sha256(canonical.encode()).hexdigest(),
            'primers_sha256': hashlib.sha256(json.dumps(primers, sort_keys=True).encode()).hexdigest()}


def tree_rss(pid):
    """Sum resident memory of the active worker and its current descendants."""
    total = 0
    try:
        text = Path('/proc/%s/status' % pid).read_text()
        rss = next(line for line in text.splitlines() if line.startswith('VmRSS:'))
        total += int(rss.split()[1]) * 1024
        children = Path('/proc/%s/task/%s/children' % (pid, pid)).read_text().split()
        total += sum(tree_rss(int(child)) for child in children)
    except (OSError, StopIteration):
        pass
    return total


def measure(args, count, mode):
    command = [sys.executable, __file__, '--worker', mode, '--count', str(count),
               '--genome', args.genome, '--db', args.db, '--threads', str(args.threads)]
    # Redirect output to files so a full pipe cannot block the worker.
    import tempfile
    with tempfile.TemporaryFile() as output, tempfile.TemporaryFile() as errors:
        process = subprocess.Popen(command, stdout=output, stderr=errors)
        peak = 0
        while process.poll() is None:
            peak = max(peak, tree_rss(process.pid))
            time.sleep(.01)
        output.seek(0); errors.seek(0)
        if process.returncode:
            raise RuntimeError(errors.read().decode(errors='replace'))
        result = json.loads(output.read())
    result['sampled_process_tree_peak_mib'] = peak / 1048576
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--genome', required=True)
    parser.add_argument('--db', required=True)
    parser.add_argument('--threads', type=int, default=2)
    parser.add_argument('--counts', default='10',
                        help='Comma-separated counts; use 10,100,500 explicitly for the expensive full panel')
    parser.add_argument('--out')
    parser.add_argument('--worker', choices=['legacy', 'batch'])
    parser.add_argument('--count', type=int)
    args = parser.parse_args()
    if args.worker:
        print(json.dumps(worker(args))); return
    report = {'database': args.db, 'genome': args.genome, 'threads': args.threads,
              'memory_method': 'Sum of /proc VmRSS for worker and descendants, sampled every 10 ms; includes BLAST, excludes supervisor. Shared resident pages may be counted more than once; this is sampled RSS, not an exact physical-memory peak.',
              'full_length_realign': True, 'results': []}
    for count in [int(value) for value in args.counts.split(',')]:
        legacy, batch = measure(args, count, 'legacy'), measure(args, count, 'batch')
        equal = legacy['evidence_sha256'] == batch['evidence_sha256']
        row = {'n_primers': count, 'legacy': legacy, 'batch': batch,
               'equal_sites_and_metadata': equal, 'speedup': legacy['seconds'] / batch['seconds']}
        report['results'].append(row)
        if args.out:
            Path(args.out).write_text(json.dumps(report, indent=2))
        print(json.dumps(row), flush=True)
        if not equal:
            raise AssertionError('Batch/legacy evidence differs')


if __name__ == '__main__':
    main()
