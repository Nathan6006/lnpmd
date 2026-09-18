"""Decide which ionizable-lipid molecules carry which protons.

Definition used throughout the project (see docs/DECISIONS.md):

    protonation level f = (protonated amines) / (all ionizable amines in the system)

With N molecules of a lipid that has s protonatable sites, the system holds
    P = f * N * s
protons. They are spread over molecules AS EVENLY AS POSSIBLE:
    every molecule gets floor(P / N) protons, and the remaining P mod N protons
    go to that many molecules chosen at random with a fixed seed.

Worked examples (N = 64):
    ECO (s=2), f=0.5 -> P=64 -> every molecule exactly +1 (no randomness)
    ECO (s=2), f=1.0 -> P=128 -> every molecule +2
    MC3 (s=1), f=0.5 -> P=32 -> 0 each, plus 32 random molecules get +1

Within one molecule, protons fill sites in `site_order` (lipids.yaml
protonation_order): a +1 ECO is protonated on site_order[0].
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np


@dataclass(frozen=True)
class MoleculeState:
    index: int                        # 0-based position of the molecule in the aggregate
    charge: int                       # number of added protons = net charge
    protonated_sites: tuple[str, ...]  # site ids that carry a proton


def n_protons(n_molecules: int, n_sites: int, fraction: float) -> int:
    """Total protons for a protonation fraction. Refuses to round."""
    if not 0.0 <= fraction <= 1.0:
        raise ValueError(f"fraction must be in [0, 1], got {fraction}")
    if n_molecules < 1 or n_sites < 1:
        raise ValueError("need at least one molecule and one site")
    exact = fraction * n_molecules * n_sites
    total = int(round(exact))
    if abs(exact - total) > 1e-9:
        raise ValueError(
            f"{fraction:.0%} of {n_molecules} molecules x {n_sites} sites = {exact} protons, "
            f"not a whole number. Change n_molecules or the protonation level.")
    return total


def assign_protonation(
    n_molecules: int,
    site_order: Sequence[str] | None,
    fraction: float,
    seed: int,
    n_sites: int | None = None,
) -> list[MoleculeState]:
    """Protonation state of every molecule.

    site_order: site ids in the order they get protonated. May be None only if
    no molecule ends up partially protonated (then order doesn't matter), in
    which case n_sites must be given.
    """
    if site_order is not None:
        n_sites = len(site_order)
    if n_sites is None:
        raise ValueError("give site_order or n_sites")

    total = n_protons(n_molecules, n_sites, fraction)
    base, extra = divmod(total, n_molecules)

    counts = np.full(n_molecules, base, dtype=int)
    if extra:
        # A dedicated generator with a fixed seed -> same molecules every time.
        rng = np.random.default_rng(seed)
        counts[rng.choice(n_molecules, size=extra, replace=False)] += 1

    states = []
    for i, k in enumerate(counts):
        k = int(k)
        if 0 < k < n_sites and site_order is None:
            raise ValueError(
                "some molecules are partially protonated, so the order in which sites are "
                "protonated matters: set protonation_order in lipids.yaml")
        sites = tuple(site_order[:k]) if site_order is not None else tuple()
        states.append(MoleculeState(i, k, sites))
    return states


def summarize(states: Sequence[MoleculeState]) -> dict[int, int]:
    """{charge: number of molecules}"""
    out: dict[int, int] = {}
    for s in states:
        out[s.charge] = out.get(s.charge, 0) + 1
    return dict(sorted(out.items()))
