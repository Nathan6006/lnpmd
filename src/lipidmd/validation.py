"""Validation: a pure POPC bilayer must reproduce published CHARMM36 numbers.

If our setup (force field files, cutoffs, thermostat, barostat, analysis
code) cannot reproduce the area per lipid, thickness and order parameters of
a plain POPC bilayer, nothing else it produces can be trusted. So:

    1. prepare  build runs/validation/<membrane>/rep1 (no aggregate)
    2. submit   run it on the cluster (same job machinery as production)
    3. analyze  membrane + order-parameter modules on the trajectory
    4. compare  against config/validation.yaml; writes validate/report.txt
                and, only if everything passes, validate/PASSED

scripts/run.py refuses to submit production runs until validate/PASSED exists.
"""

from __future__ import annotations

import copy
import datetime as dt
from pathlib import Path

import pandas as pd

from .config import Config, ConfigError


def validation_config(cfg: Config) -> Config:
    """Copy of cfg with the validation run's production length."""
    v = copy.copy(cfg)
    v.data = copy.deepcopy(cfg.data)
    v.data["md"]["production"]["length_ns"] = cfg.require("validation.production_ns")
    return v


def run_dir(cfg: Config, smoke: bool = False) -> Path:
    base = cfg.systems["smoke_runs_dir"] if smoke else cfg.systems["runs_dir"]
    return cfg.path(base) / "validation" / cfg.validation["membrane"] / "rep1"


def results_dir(cfg: Config) -> Path:
    return cfg.path(cfg.systems["results_dir"]) / "validation" / cfg.validation["membrane"] / "rep1"


def marker(cfg: Config) -> Path:
    return cfg.path("validate") / "PASSED"


def compare(cfg: Config) -> tuple[bool, str]:
    v = cfg.validation
    discard_ps = float(cfg.require("validation.discard_ns")) * 1000
    res = results_dir(cfg)
    memb = pd.read_csv(res / "membrane.csv")
    memb = memb[memb.time_ps >= discard_ps]
    order = pd.read_csv(res / "order.csv")
    order = order[(order.time_ps >= discard_ps) & (order.resname == "POPC")]
    if memb.empty:
        raise ValueError(f"no membrane data after discard_ns = {v['discard_ns']}")

    sim = {
        "area_per_lipid_nm2": memb.loc[memb.observable == "area_per_lipid_nm2", "value"].mean(),
        "thickness_dpp_nm": memb.loc[memb.observable == "thickness_dpp_nm", "value"].mean(),
    }
    scd = order.groupby(["chain", "carbon"])["scd"].mean()

    lines = [f"Validation report ({dt.datetime.now():%Y-%m-%d %H:%M}) -- membrane {v['membrane']}",
             f"discarding first {v['discard_ns']} ns", ""]
    ok = True
    problems = []
    for key in ("area_per_lipid_nm2", "thickness_dpp_nm"):
        ref = v["references"][key]
        if ref["value"] is None or ref["tolerance"] is None:
            problems.append(f"references.{key}: value/tolerance not filled in")
            lines.append(f"  {key:22s} sim {sim[key]:.4f}   ref: NOT SET")
            ok = False
            continue
        diff = sim[key] - float(ref["value"])
        passed = abs(diff) <= float(ref["tolerance"])
        ok &= passed
        lines.append(f"  {key:22s} sim {sim[key]:.4f}  ref {ref['value']}  diff {diff:+.4f}  "
                     f"tol {ref['tolerance']}  {'PASS' if passed else 'FAIL'}")
    for key, chain in (("scd_sn1", "sn1"), ("scd_sn2", "sn2")):
        ref = v["references"][key]
        if not ref.get("values") or ref.get("tolerance") is None:
            problems.append(f"references.{key}: values/tolerance not filled in")
            lines.append(f"  {key:22s} ref: NOT SET")
            ok = False
            continue
        for carbon, rv in ref["values"].items():
            if (chain, carbon) not in scd.index:
                lines.append(f"  {key} {carbon}: not computed (check tails in analysis.yaml)")
                ok = False
                continue
            s = scd[(chain, carbon)]
            passed = abs(s - float(rv)) <= float(ref["tolerance"])
            ok &= passed
            lines.append(f"  {key} {carbon:5s} sim {s:+.4f}  ref {float(rv):+.4f}  {'PASS' if passed else 'FAIL'}")
    if problems:
        lines += ["", "Cannot pass until config/validation.yaml is complete:"] + [f"  - {p}" for p in problems]
    lines += ["", "RESULT: " + ("PASS" if ok else "FAIL")]
    text = "\n".join(lines) + "\n"

    out = cfg.path("validate")
    out.mkdir(exist_ok=True)
    (out / "report.txt").write_text(text)
    if ok:
        marker(cfg).write_text(text)
    else:
        marker(cfg).unlink(missing_ok=True)
    return ok, text


def require_passed(cfg: Config) -> None:
    if not marker(cfg).exists():
        raise ConfigError("validation has not passed (validate/PASSED missing). Run the validation "
                          "first -- README step 5 -- or pass --skip-validation-check if you really mean it.")
