#!/usr/bin/env python3
"""Summary table: rows = 9 systems, columns = observables, mean ± SD over replicas.

    python scripts/report.py

Reads results/<system>/rep<N>/*.csv; writes results/summary_{long,table}.csv
and results/summary_table.md.
"""
import _setup  # noqa: F401
import argparse
import sys

from lipidmd.config import Config, ConfigError
from lipidmd.report import write_report


def main() -> int:
    argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter).parse_args()
    try:
        out = write_report(Config())
    except ConfigError as e:
        print(e)
        return 1
    print((out / "summary_table.md").read_text())
    print(f"written to {out}/summary_table.csv (+ summary_long.csv, summary_table.md)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
