"""Compare bounded synthetic FASTA reads with and without a reused handle."""
import argparse
import json
from pathlib import Path
import platform
import tempfile
import time

from primerblast_oss.genome import Genome


def run(reads=500, work_dir=None):
    if reads < 1:
        raise ValueError('reads must be positive')
    if work_dir:
        Path(work_dir).mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='synthetic-fasta-', dir=work_dir) as folder:
        path = Path(folder) / 'synthetic.fa'
        sequence = 'ACGTTGCAAGTCCGATCGTA' * 200
        header = '>synthetic\n'
        path.write_bytes((header + '\n'.join(sequence[i:i + 80] for i in range(0, len(sequence), 80)) + '\n').encode('ascii'))
        Path(str(path) + '.fai').write_text('synthetic\t{}\t{}\t80\t81\n'.format(len(sequence), len(header)), encoding='ascii')
        genome = Genome(str(path))
        positions = [1 + i * 37 % (len(sequence) - 26) for i in range(reads)]
        started = time.perf_counter()
        separate = [genome.fetch('synthetic', start, start + 25) for start in positions]
        separate_seconds = time.perf_counter() - started
        started = time.perf_counter()
        with genome.read_session():
            reused = [genome.fetch('synthetic', start, start + 25) for start in positions]
        reused_seconds = time.perf_counter() - started
        return {'benchmark': 'synthetic_read_session', 'python': platform.python_version(),
                'reads': reads, 'window_bases': 26, 'reference_bases': len(sequence),
                'identical_bases': separate == reused,
                'separate_open_seconds': round(separate_seconds, 6),
                'read_session_seconds': round(reused_seconds, 6)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reads', type=int, default=500)
    parser.add_argument('--work-dir')
    parser.add_argument('--json-out')
    args = parser.parse_args()
    result = run(args.reads, args.work_dir)
    output = json.dumps(result, indent=2, sort_keys=True) + '\n'
    if args.json_out:
        Path(args.json_out).write_text(output, encoding='utf-8')
    print(output, end='')
    return 0 if result['identical_bases'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
