"""Convergence diagnostics for scalar time series.

MD frames are strongly correlated in time, so the naive standard error
(std / sqrt(n_frames)) is far too small. Two diagnostics:

running average   mean of everything up to time t. If it is still drifting
                  at the end of the run, the run is not long enough (or the
                  system is still equilibrating -- see discard_ns).

block averaging   (Flyvbjerg & Petersen, J. Chem. Phys. 91, 461, 1989)
                  Split the series into blocks of 1, 2, 4, 8, ... frames and
                  compute the standard error of the block means. While blocks
                  are shorter than the correlation time the estimate grows;
                  once blocks are longer, it levels off at the true standard
                  error. A plateau = trustworthy error bar. No plateau = the
                  run is too short to estimate its own error.

This module reads the per-frame CSVs the other modules wrote, so run it last.
Outputs: convergence_running.csv, convergence_blocking.csv.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .common import AnalysisContext


def running_average(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, float)
    return np.cumsum(x) / np.arange(1, len(x) + 1)


def block_standard_errors(x: np.ndarray, min_blocks: int = 4) -> pd.DataFrame:
    """Standard error of the mean vs block size (in frames), by repeated halving."""
    x = np.asarray(x, float)
    rows = []
    size = 1
    while len(x) // size >= min_blocks:
        n_blocks = len(x) // size
        means = x[: n_blocks * size].reshape(n_blocks, size).mean(axis=1)
        se = means.std(ddof=1) / np.sqrt(n_blocks)
        rows.append((size, n_blocks, se))
        size *= 2
    return pd.DataFrame(rows, columns=["block_frames", "n_blocks", "sem"])


def scalar_series(ctx: AnalysisContext) -> dict[str, pd.DataFrame]:
    """Collect (time_ps, value) series for every per-frame scalar observable."""
    out: dict[str, pd.DataFrame] = {}
    d = ctx.out_dir
    for name in ("membrane", "sasa", "water_core"):
        p = d / f"{name}.csv"
        if p.exists():
            df = pd.read_csv(p)
            for obs, g in df.groupby("observable"):
                out[obs] = g[["time_ps", "value"]].reset_index(drop=True)
    p = d / "insertion.csv"
    if p.exists():
        df = pd.read_csv(p)
        out["insertion_com_dz_mean_nm"] = df.groupby("time_ps")["com_dz_nm"].mean().rename("value").reset_index()
    p = d / "contacts.csv"
    if p.exists():
        df = pd.read_csv(p)
        if len(df):
            out["n_amines_in_contact"] = (df.groupby("time_ps")["n_amines_in_contact"].sum()
                                          .rename("value").reset_index())
    p = d / "order.csv"
    if p.exists():
        df = pd.read_csv(p)
        for (res, chain), g in df.groupby(["resname", "chain"]):
            out[f"scd_mean_{res}_{chain}"] = g.groupby("time_ps")["scd"].mean().rename("value").reset_index()
    return out


def compute(ctx: AnalysisContext) -> dict[str, pd.DataFrame]:
    run_rows, blk_rows = [], []
    for obs, s in scalar_series(ctx).items():
        ra = running_average(s["value"].to_numpy())
        run_rows.append(pd.DataFrame({"observable": obs, "time_ps": s["time_ps"], "value": s["value"],
                                      "running_mean": ra}))
        b = block_standard_errors(s["value"].to_numpy())
        b.insert(0, "observable", obs)
        blk_rows.append(b)
    cols_r = ["observable", "time_ps", "value", "running_mean"]
    cols_b = ["observable", "block_frames", "n_blocks", "sem"]
    return {
        "convergence_running": pd.concat(run_rows, ignore_index=True) if run_rows else pd.DataFrame(columns=cols_r),
        "convergence_blocking": pd.concat(blk_rows, ignore_index=True) if blk_rows else pd.DataFrame(columns=cols_b),
    }
