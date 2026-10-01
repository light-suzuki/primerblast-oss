"""Launch a synthetic-only GUI for reproducible documentation screenshots.

Run from a source checkout with BLAST+ and Primer3 on PATH. Existing reference
catalogs and databases are not loaded. Generated files stay in .local by default.
"""
import argparse
import json
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from benchmarks.continuous_benchmark import build_fixture
from primerblast_oss.webapp import server


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--port', type=int, default=8880)
    parser.add_argument('--workdir', default=str(ROOT / '.local/readme-demo'))
    args = parser.parse_args()
    demo = Path(args.workdir).resolve()
    demo.mkdir(parents=True, exist_ok=True)
    os.chdir(demo)
    databases = Path('databases')
    databases.mkdir(exist_ok=True)
    os.environ['BLAST_USAGE_REPORT'] = 'false'
    build_fixture(databases)
    (databases / 'synthetic.gff3').write_text(
        '##gff-version 3\nchr_target\tfixture\tgene\t10\t80\t.\t-\t.\tID=g1;Name=ExampleGene\n',
        encoding='utf-8')
    Path('references.json').write_text(json.dumps([{
        'name': 'Synthetic PCR demo', 'genome': 'databases/synthetic.fa',
        'gff3': 'databases/synthetic.gff3', 'database': 'databases/synthetic_db',
    }]), encoding='utf-8')
    os.environ['PRIMERBLAST_REFERENCES'] = 'references.json'
    server.DEFAULT_DB_DIRS = [databases]
    httpd = server.serve('127.0.0.1', args.port)
    print('Synthetic-only GUI: http://localhost:{}/'.format(args.port), flush=True)
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        httpd.server_close()


if __name__ == '__main__':
    main()
