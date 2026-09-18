#!/usr/bin/env python3
"""Load every config file, check it is consistent, and list open decisions.

    python scripts/check_config.py

Every `null` in config/ is either a fact you still need to fill in (TODO) or
a scientific decision nobody has made yet (ASK). Nothing silently defaults.
"""
import _setup  # noqa: F401
import argparse
import sys
from collections import defaultdict

from lipidmd.config import Config, ConfigError


def main() -> int:
    argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter).parse_args()
    try:
        cfg = Config()
    except ConfigError as e:
        print(f"CONFIG ERROR: {e}")
        return 2
    print(f"Config OK: {len(cfg.systems_list())} systems x {cfg.systems['replicas']} replicas "
          f"= {len(cfg.runs())} runs; target membrane = {cfg.target_membrane()}")
    for s in cfg.systems_list():
        print(f"  {s.name}")
    try:
        counts = cfg.membrane_counts()
        print(f"Per-leaflet composition: {counts}")
    except ConfigError as e:
        print(f"Membrane composition: {e}")

    open_items = cfg.unresolved()
    if not open_items:
        print("\nNo unresolved values.")
        return 0
    by_file = defaultdict(list)
    for p in open_items:
        f, _, rest = p.partition(".")
        by_file[f].append(rest)
    print(f"\n{len(open_items)} values are still null (see the comment beside each in the file,")
    print("and docs/DECISIONS.md 'Open questions'):")
    for f, items in by_file.items():
        print(f"\n  config/{f}.yaml")
        for it in items:
            print(f"    - {it}")
    return 1


if __name__ == "__main__":
    sys.exit(main())
