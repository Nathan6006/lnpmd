"""Contacts between ionizable-lipid amine N atoms and POPS anionic oxygens.

A contact between one amine N and one POPS molecule exists in a frame if any
of that POPS's phosphate/carboxylate oxygens (analysis.yaml
contacts.pops_oxygens) is within `cutoff_nm` of the N (periodic boundaries
respected).

Outputs
    contacts.csv         per frame: how many amine N are in contact, split by
                         the charge of the molecule they belong to, and the
                         number of N-POPS pairs
    contact_events.csv   one row per contact event (N atom, POPS residue,
                         start, duration). Residence time = duration.
                         Breaks of <= gap_tolerance_frames are bridged.
                         Events still open at the end are flagged `censored`
                         (their true duration is longer than recorded).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .common import A_PER_NM, AnalysisContext


class EventTracker:
    """Turns per-frame sets of contacting pairs into contact events.

    Kept separate from MDAnalysis so it can be unit tested with plain sets.
    """

    def __init__(self, gap_tolerance: int):
        self.gap = gap_tolerance
        self.open: dict[tuple, list] = {}   # pair -> [start_frame, last_seen_frame]
        self.events: list[tuple] = []       # (pair, start_frame, last_seen_frame, censored)

    def update(self, frame: int, pairs: set[tuple]) -> None:
        for p in pairs:
            if p in self.open:
                self.open[p][1] = frame
            else:
                self.open[p] = [frame, frame]
        for p in list(self.open):
            start, last = self.open[p]
            if frame - last > self.gap:
                self.events.append((p, start, last, False))
                del self.open[p]

    def close(self) -> None:
        for p, (start, last) in self.open.items():
            self.events.append((p, start, last, True))
        self.open.clear()


def compute(ctx: AnalysisContext) -> dict[str, pd.DataFrame]:
    from MDAnalysis.lib.distances import capped_distance

    c = ctx.cfg.analysis["contacts"]
    if c.get("cutoff_nm") is None or c.get("gap_tolerance_frames") is None:
        raise ValueError("analysis.yaml contacts.cutoff_nm / gap_tolerance_frames are null (open decisions)")
    cutoff = float(c["cutoff_nm"]) * A_PER_NM
    u = ctx.universe
    amines = ctx.amine_atoms()
    oxy = u.select_atoms(f"resname POPS and name {' '.join(c['pops_oxygens'])}")
    if oxy.n_atoms == 0:
        # e.g. pure POPC control: no anionic lipid, so no contacts by definition
        return {"contacts": pd.DataFrame(columns=["time_ps", "charge", "n_amines_in_contact", "n_amines", "n_pairs"]),
                "contact_events": pd.DataFrame(columns=["amine_atom", "pops_residue_index", "start_ps", "duration_ps", "censored"])}
    q_of = ctx.charge_of_resname()
    amine_charge = np.array([q_of[r] for r in amines.resnames])
    ox_resid = oxy.resindices  # 0-based residue index, unique even in huge systems

    idx_of = {a.index: n for n, a in enumerate(amines)}
    tracker = EventTracker(int(c["gap_tolerance_frames"]))
    rows, times = [], []
    for k, ts in enumerate(ctx.frames()):
        times.append(ts.time)
        pairs_idx = capped_distance(amines.positions, oxy.positions, cutoff, box=ts.dimensions,
                                    return_distances=False)
        pairs = {(int(amines[i].index), int(ox_resid[j])) for i, j in pairs_idx}
        tracker.update(k, pairs)
        in_contact = np.zeros(amines.n_atoms, bool)
        for a, _ in pairs:
            in_contact[idx_of[a]] = True
        for q in np.unique(amine_charge):
            m = amine_charge == q
            rows.append((ts.time, int(q), int(in_contact[m].sum()), int(m.sum()),
                         sum(1 for a, _ in pairs if amine_charge[idx_of[a]] == q)))
    tracker.close()

    dt = (times[1] - times[0]) if len(times) > 1 else 0.0
    ev = pd.DataFrame([(p[0], p[1], times[s], (e - s + 1) * dt, cen) for p, s, e, cen in tracker.events],
                      columns=["amine_atom", "pops_residue_index", "start_ps", "duration_ps", "censored"])
    ts_df = pd.DataFrame(rows, columns=["time_ps", "charge", "n_amines_in_contact", "n_amines", "n_pairs"])
    return {"contacts": ts_df, "contact_events": ev}
