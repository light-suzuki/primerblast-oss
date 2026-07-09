"""Launch the primerblast-oss web GUI.

    python -m primerblast_oss.webapp [--host H] [--port P] [--no-browser]

Binds to 127.0.0.1 by default (local single-user use). BLAST/primer3 run on
whatever machine hosts this process, so start it where the databases live
(e.g. inside WSL).
"""
from __future__ import annotations

import argparse
import sys
import threading
import webbrowser

from .server import serve, health


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="primerblast-oss webapp")
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--port", type=int, default=8799)
    ap.add_argument("--no-browser", action="store_true",
                    help="do not try to open a browser tab")
    args = ap.parse_args(argv)

    hc = health()
    if not hc["ok"]:
        missing = [k for k, v in hc["tools"].items() if not v]
        print(f"[warning] required tools missing: {', '.join(missing)}. "
              "Design/specificity will fail until primer3_core and blastn are on PATH.",
              file=sys.stderr)

    httpd = serve(args.host, args.port)
    url = f"http://{args.host}:{args.port}/"
    print(f"primerblast-oss GUI  ->  {url}")
    print("Press Ctrl+C to stop.")

    if not args.no_browser:
        threading.Timer(0.6, lambda: webbrowser.open(url)).start()

    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\nshutting down...")
    finally:
        httpd.server_close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
