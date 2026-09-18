"""Collect per-replica results into one table: rows = systems, columns =
observables, mean +/- SD across replicas.

Per replica, each observable is first reduced to ONE number (its time
average after discarding the first `discard_ns` of production). The table
then reports the mean and standard deviation of those replica numbers. The
SD across independent replicas is the honest error bar here: it includes
everything that differs between runs, which frame-to-frame scatter does not.

Outputs (in results/):
    summary_long.csv    system, lipid, protonation, replica, observable, value
    summary_table.csv   one row per system, <obs>_mean / <obs>_sd / <obs>_n
    summary_table.md    the same, formatted "mean ± SD", for pasting
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .config import Config


def _after(df: pd.DataFrame, t0_ps: float) -> pd.DataFrame:
    return df[df["time_ps"] >= t0_ps]


def replica_scalars(res_dir: Path, discard_ns: float) -> dict[str, float]:
    """One number per observable for one replica's results directory."""
    t0 = discard_ns * 1000.0
    out: dict[str, float] = {}

    def read(name: str) -> pd.DataFrame | None:
        p = res_dir / f"{name}.csv"
        return pd.read_csv(p) if p.exists() else None

    if (df := read("membrane")) is not None:
        df = _after(df, t0)
        for obs in ("area_per_lipid_nm2", "thickness_dpp_nm"):
            out[obs] = df.loc[df.observable == obs, "value"].mean()
    if (df := read("order")) is not None:
        df = _after(df, t0)
        for (res, chain), g in df.groupby(["resname", "chain"]):
            out[f"scd_{res}_{chain}_mean"] = g["scd"].mean()
    if (df := read("insertion")) is not None:
        df = _after(df, t0)
        out["insertion_com_dz_nm"] = df["com_dz_nm"].mean()
        out["insertion_amine_dz_nm"] = df["amine_dz_nm"].mean()
    if (df := read("contacts")) is not None and len(df):
        df = _after(df, t0)
        out["amines_in_contact"] = df.groupby("time_ps")["n_amines_in_contact"].sum().mean()
    if (df := read("contact_events")) is not None and len(df):
        df = df[df["start_ps"] >= t0]
        # censored events (still in contact at the end) are included at their
        # observed length, which biases this mean LOW; see DECISIONS.md.
        out["contact_residence_ps"] = df["duration_ps"].mean()
    if (df := read("sasa")) is not None:
        df = _after(df, t0)
        out["sasa_hydrophobic_nm2"] = df.loc[df.observable == "sasa_hydrophobic_nm2", "value"].mean()
    if (df := read("water_core")) is not None:
        out["waters_in_core"] = _after(df, t0)["value"].mean()
    if (df := read("water_permeation")) is not None:
        out["permeation_events"] = float((df["time_ps"] >= t0).sum())
    return out


def build_report(cfg: Config, results_dir: Path | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    discard = float(cfg.require("analysis.discard_ns"))
    results_dir = results_dir or cfg.path(cfg.systems["results_dir"])
    rows = []
    for run in cfg.runs():
        d = results_dir / run.system.name / f"rep{run.replica}"
        if not d.exists():
            continue
        for obs, v in replica_scalars(d, discard).items():
            rows.append((run.system.name, run.system.lipid, run.system.fraction, run.replica, obs, v))
    long = pd.DataFrame(rows, columns=["system", "lipid", "protonation", "replica", "observable", "value"])
    if long.empty:
        return long, pd.DataFrame()
    agg = (long.groupby(["system", "lipid", "protonation", "observable"])["value"]
           .agg(mean="mean", sd=lambda x: x.std(ddof=1) if len(x) > 1 else np.nan, n="count")
           .reset_index())
    wide = agg.pivot_table(index=["system", "lipid", "protonation"], columns="observable",
                           values=["mean", "sd", "n"])
    wide.columns = [f"{obs}_{stat}" for stat, obs in wide.columns]
    wide = wide[sorted(wide.columns)].reset_index()
    order = {s.name: i for i, s in enumerate(cfg.systems_list())}
    wide = wide.sort_values("system", key=lambda s: s.map(order)).reset_index(drop=True)
    return long, wide


def to_markdown(wide: pd.DataFrame) -> str:
    obs = sorted({c.rsplit("_", 1)[0] for c in wide.columns if c.endswith("_mean")})
    head = "| system | " + " | ".join(obs) + " |"
    sep = "|---" * (len(obs) + 1) + "|"
    lines = [head, sep]
    for _, r in wide.iterrows():
        cells = []
        for o in obs:
            m, s, n = r.get(f"{o}_mean"), r.get(f"{o}_sd"), r.get(f"{o}_n")
            if pd.isna(m):
                cells.append("—")
            else:
                sd = "n/a" if pd.isna(s) else f"{s:.3g}"
                cells.append(f"{m:.4g} ± {sd} (n={int(n)})")
        lines.append(f"| {r['system']} | " + " | ".join(cells) + " |")
    return "\n".join(lines) + "\n"


def write_report(cfg: Config) -> Path:
    out = cfg.path(cfg.systems["results_dir"])
    out.mkdir(parents=True, exist_ok=True)
    long, wide = build_report(cfg, out)
    long.to_csv(out / "summary_long.csv", index=False)
    wide.to_csv(out / "summary_table.csv", index=False)
    (out / "summary_table.md").write_text(to_markdown(wide) if not wide.empty else "no results yet\n")
    return out
