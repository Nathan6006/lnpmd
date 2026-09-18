"""One module per observable. Each exposes

    compute(ctx: AnalysisContext) -> dict[str, pandas.DataFrame]

returning one or more tidy tables (one row per observation) that are written
to results/<system>/rep<N>/<table name>.csv.

Run them with scripts/analyze.py.
"""

OBSERVABLES = ["membrane", "order", "insertion", "contacts", "sasa", "water"]
