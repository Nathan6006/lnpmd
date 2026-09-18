"""CGenFF parameter ingestion and checking.

Background (why penalties matter)
---------------------------------
CGenFF assigns parameters to a new molecule BY ANALOGY to molecules it was
fitted on. Each assigned parameter and each partial charge gets a "penalty"
saying how far the analogy was stretched. The CGenFF authors' guidance:
    penalty < 10     analogy is fair, generally usable as is
    10 <= p <= 50    usable, but do some basic validation
    p > 50           poor analogy; the parameter should be optimized/validated
                     before it can be trusted
Penalties near the ionizable nitrogens matter most for this project, because
the N and its neighbours (their charges, and the torsions that set where the
headgroup points) are exactly the chemistry whose protonation we compare.

What a .str file looks like (abridged)::

    RESI MC3H    1.000 ! param penalty=  24.500 ; charge penalty=  31.200
    ATOM N1     NG3P1  -0.310 !   12.345      <- per-atom charge penalty
    BOND N1   C1   N1   C2  ...
    read param card flex append
    BONDS
    CG321  NG3P1   200.00   1.4800 ! MC3H , from CG321 NG3P0, penalty= 4.0
    DIHEDRALS
    CG321  CG321  CG321  NG3P1   0.1000  3   0.00 ! MC3H , from ..., penalty= 61.5
"""

from __future__ import annotations

import re
import shutil
import subprocess
from collections import Counter, deque
from dataclasses import dataclass, field
from pathlib import Path

REVIEW_THRESHOLD = 10.0
ATTENTION_THRESHOLD = 50.0

_PARAM_SECTIONS = {"BONDS": ("bond", 2), "ANGLES": ("angle", 3), "DIHEDRALS": ("dihedral", 4),
                   "IMPROPERS": ("improper", 4), "IMPROPER": ("improper", 4)}
_PENALTY = re.compile(r"penalty\s*=\s*([-\d.]+)", re.IGNORECASE)


@dataclass
class StrMolecule:
    """What we need from a CGenFF .str file."""

    resname: str
    net_charge: float
    param_penalty: float | None
    charge_penalty: float | None
    atom_names: list[str] = field(default_factory=list)
    atom_types: dict[str, str] = field(default_factory=dict)
    atom_charges: dict[str, float] = field(default_factory=dict)
    atom_charge_penalties: dict[str, float] = field(default_factory=dict)
    bonds: list[tuple[str, str]] = field(default_factory=list)
    impropers: list[tuple[str, str, str, str]] = field(default_factory=list)
    parameters: list["Parameter"] = field(default_factory=list)


@dataclass
class Parameter:
    kind: str                 # bond | angle | dihedral | improper
    types: tuple[str, ...]
    penalty: float
    line: str


@dataclass
class PenaltyItem:
    kind: str                 # bond | angle | dihedral | improper | charge
    label: str                # atom types (parameters) or atom name (charges)
    penalty: float
    instances: list[tuple[str, ...]]   # atom-name tuples where this applies
    near_amine: bool
    line: str = ""

    @property
    def severity(self) -> str:
        return "ATTENTION" if self.penalty > ATTENTION_THRESHOLD else "review"


def parse_str(path: str | Path) -> StrMolecule:
    mol: StrMolecule | None = None
    mode = None           # "rtf" | "param"
    section = None
    for raw in Path(path).read_text().splitlines():
        line = raw.strip()
        up = line.upper()
        if not line or line.startswith("*"):
            continue
        if up.startswith("READ RTF"):
            mode = "rtf"
            continue
        if up.startswith("READ PARA"):
            mode, section = "param", None
            continue
        if up.startswith("END") or up.startswith("RETURN"):
            mode = None if up.startswith("RETURN") else mode
            section = None
            continue

        body, _, comment = line.partition("!")
        f = body.split()
        if not f:
            continue
        key = f[0].upper()

        if mode == "rtf":
            if key == "RESI":
                pp = re.search(r"param penalty\s*=\s*([-\d.]+)", comment)
                cp = re.search(r"charge penalty\s*=\s*([-\d.]+)", comment)
                mol = StrMolecule(f[1], float(f[2]),
                                  float(pp.group(1)) if pp else None,
                                  float(cp.group(1)) if cp else None)
            elif key == "ATOM" and mol is not None:
                name, atype, q = f[1], f[2], float(f[3])
                mol.atom_names.append(name)
                mol.atom_types[name] = atype
                mol.atom_charges[name] = q
                if comment.strip():
                    try:
                        mol.atom_charge_penalties[name] = float(comment.split()[0])
                    except ValueError:
                        pass
            elif key in ("BOND", "DOUBLE", "TRIPLE") and mol is not None:
                names = f[1:]
                mol.bonds += [(names[i], names[i + 1]) for i in range(0, len(names) - 1, 2)]
            elif key in ("IMPR", "IMPH") and mol is not None:
                names = f[1:]
                mol.impropers += [tuple(names[i:i + 4]) for i in range(0, len(names) - 3, 4)]
        elif mode == "param":
            if key in _PARAM_SECTIONS and len(f) == 1:
                section = key
                continue
            if key in ("NONBONDED", "CMAP", "NBFIX", "HBOND"):
                section = None
                continue
            if section is None or mol is None:
                continue
            kind, n = _PARAM_SECTIONS[section]
            m = _PENALTY.search(comment)
            mol.parameters.append(Parameter(kind, tuple(f[:n]), float(m.group(1)) if m else 0.0, line))
    if mol is None:
        raise ValueError(f"{path}: no RESI record found -- is this a CGenFF .str file?")
    return mol


# --------------------------------------------------------------------------
# Graph helpers: which atoms are "near" the ionizable nitrogens?
# --------------------------------------------------------------------------
def _adjacency(mol: StrMolecule) -> dict[str, set[str]]:
    adj: dict[str, set[str]] = {n: set() for n in mol.atom_names}
    for a, b in mol.bonds:
        adj.setdefault(a, set()).add(b)
        adj.setdefault(b, set()).add(a)
    return adj


def atoms_within(mol: StrMolecule, centers: list[str], max_bonds: int) -> set[str]:
    """Atoms at most `max_bonds` bonds away from any center atom (BFS)."""
    adj = _adjacency(mol)
    seen = {c: 0 for c in centers if c in adj}
    queue = deque(seen)
    while queue:
        a = queue.popleft()
        if seen[a] == max_bonds:
            continue
        for b in adj[a]:
            if b not in seen:
                seen[b] = seen[a] + 1
                queue.append(b)
    return set(seen)


def _instances(mol: StrMolecule, kind: str) -> list[tuple[str, ...]]:
    """Every bond/angle/dihedral/improper in the molecule as atom-name tuples."""
    adj = _adjacency(mol)
    if kind == "bond":
        return [tuple(b) for b in mol.bonds]
    if kind == "angle":
        return [(a, b, c) for b in adj for a in adj[b] for c in adj[b] if a < c]
    if kind == "dihedral":
        out = []
        for b, c in mol.bonds:
            for a in adj[b] - {c}:
                for d in adj[c] - {b, a}:
                    out.append((a, b, c, d))
        return out
    if kind == "improper":
        return list(mol.impropers)
    raise ValueError(kind)


def _matches(types: tuple[str, ...], inst: tuple[str, ...], mol: StrMolecule) -> bool:
    t = tuple(mol.atom_types.get(a, "?") for a in inst)
    return t == types or t[::-1] == types


def check_penalties(mol: StrMolecule, amine_atoms: list[str], threshold: float = REVIEW_THRESHOLD,
                    near_bonds: int = 3) -> list[PenaltyItem]:
    """Every parameter or charge with penalty > threshold, most severe first.

    amine_atoms: names of the protonatable nitrogens. An item is flagged
    `near_amine` if any atom it applies to is within `near_bonds` bonds of one
    of them.
    """
    near = atoms_within(mol, amine_atoms, near_bonds) if amine_atoms else set()
    items: list[PenaltyItem] = []
    cache: dict[str, list[tuple[str, ...]]] = {}
    for p in mol.parameters:
        if p.penalty <= threshold:
            continue
        pool = cache.setdefault(p.kind, _instances(mol, p.kind))
        inst = [i for i in pool if _matches(p.types, i, mol)]
        items.append(PenaltyItem(p.kind, " ".join(p.types), p.penalty, inst,
                                  any(a in near for i in inst for a in i), p.line))
    for name, pen in mol.atom_charge_penalties.items():
        if pen > threshold:
            items.append(PenaltyItem("charge", f"{name} ({mol.atom_types[name]})", pen,
                                     [(name,)], name in near))
    items.sort(key=lambda it: (-it.penalty, it.kind))
    return items


def format_penalty_report(mol: StrMolecule, items: list[PenaltyItem], amine_atoms: list[str],
                          threshold: float = REVIEW_THRESHOLD, near_bonds: int = 3) -> str:
    out = [f"CGenFF penalty report for {mol.resname} (net charge {mol.net_charge:+.3f})",
           f"  overall: param penalty = {mol.param_penalty}, charge penalty = {mol.charge_penalty}",
           f"  showing penalties > {threshold:g}; > {ATTENTION_THRESHOLD:g} is flagged ATTENTION"]
    if amine_atoms:
        out.append(f"  '*' = within {near_bonds} bonds of an ionizable N ({', '.join(amine_atoms)})")
    else:
        out.append("  WARNING: protonatable atom names unknown (lipids.yaml), so proximity to the "
                   "amine is NOT checked")
    missing = [a for a in amine_atoms if a not in mol.atom_types]
    if missing:
        out.append(f"  WARNING: amine atom(s) {missing} not found in this .str -- check lipids.yaml")
    if abs(mol.net_charge - round(mol.net_charge)) > 1e-3:
        out.append(f"  WARNING: net charge {mol.net_charge} is not an integer")
    if not items:
        out.append("  nothing above threshold.")
        return "\n".join(out)
    out.append("")
    out.append(f"  {'':1s} {'severity':9s} {'kind':9s} {'penalty':>8s}  types / atom       applies to")
    for it in items:
        where = "; ".join("-".join(i) for i in it.instances[:3])
        if len(it.instances) > 3:
            where += f"; ... ({len(it.instances)} total)"
        if not it.instances:
            where = "(no matching atoms found)"
        out.append(f"  {'*' if it.near_amine else ' '} {it.severity:9s} {it.kind:9s} "
                   f"{it.penalty:8.1f}  {it.label:18s} {where}")
    n_att = sum(it.penalty > ATTENTION_THRESHOLD for it in items)
    n_near = sum(it.near_amine for it in items)
    out += ["", f"  {len(items)} above threshold, {n_att} need manual attention, "
                f"{n_near} near an ionizable N."]
    if any(it.near_amine and it.penalty > ATTENTION_THRESHOLD for it in items):
        out.append("  >>> High-penalty parameters touch the ionizable headgroup. Do not run "
                   "production with these until they are validated (see DECISIONS.md).")
    return "\n".join(out)


# --------------------------------------------------------------------------
# Structure sanity checks
# --------------------------------------------------------------------------
def mol2_formula(path: str | Path) -> str:
    """Hill-order formula (C, H, then alphabetical) from a .mol2's atom block."""
    counts: Counter[str] = Counter()
    in_atoms = False
    for raw in Path(path).read_text().splitlines():
        if raw.startswith("@<TRIPOS>"):
            in_atoms = raw.strip() == "@<TRIPOS>ATOM"
            continue
        f = raw.split()
        if in_atoms and len(f) >= 6:
            # column 6 is the SYBYL type, e.g. "C.3", "N.4", "H", "O.co2"
            counts[f[5].split(".")[0].capitalize()] += 1
    return hill_formula(counts)


def hill_formula(counts: dict[str, int]) -> str:
    order = [e for e in ("C", "H") if e in counts] + sorted(e for e in counts if e not in ("C", "H"))
    return "".join(f"{e}{counts[e] if counts[e] != 1 else ''}" for e in order)


def parse_formula(formula: str) -> Counter[str]:
    out: Counter[str] = Counter()
    for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
        out[el] += int(n) if n else 1
    return out


def expected_formula(neutral: str, n_added_h: int) -> str:
    c = parse_formula(neutral)
    c["H"] += n_added_h
    return hill_formula(c)


# --------------------------------------------------------------------------
# Paths + conversion to GROMACS format
# --------------------------------------------------------------------------
def cgenff_files(cfg, lipid: str, charge: int) -> dict[str, Path]:
    resname = cfg.lipid(lipid)["states"][charge]["resname"]
    cg = cfg.path(cfg.paths["cgenff_dir"]) / lipid
    pr = cfg.path(cfg.paths["params_dir"]) / lipid
    return {"mol2": cg / f"{resname}.mol2", "str": cg / f"{resname}.str",
            "itp": pr / f"{resname}.itp", "prm": pr / f"{resname}.prm", "resname": resname}


def find_param_file(path: Path) -> Path:
    """cgenff_charmm2gmx names outputs in lowercase; accept either spelling."""
    if path.exists():
        return path
    low = path.with_name(path.name.lower())
    if low.exists():
        return low
    raise FileNotFoundError(f"{path} (or {low.name}) not found -- run scripts/prepare_params.py")


def convert_to_gromacs(cfg, lipid: str, charge: int) -> list[str]:
    """Run the MacKerell-lab cgenff_charmm2gmx script for one charge state.

    Returns log lines. The script writes <resname>.itp/.prm/_ini.pdb into the
    current directory, so we run it inside the params folder.
    """
    files = cgenff_files(cfg, lipid, charge)
    script, ff = cfg.paths.get("cgenff_charmm2gmx"), cfg.paths.get("charmm36_ff_dir")
    if not script or not ff:
        raise RuntimeError("set cgenff_charmm2gmx and charmm36_ff_dir in config/paths.yaml "
                           "(MANUAL_STEPS.md step 3)")
    for k in ("mol2", "str"):
        if not files[k].exists():
            raise FileNotFoundError(files[k])
    out_dir = files["itp"].parent
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = ["python3", str(cfg.path(script)), files["resname"], str(files["mol2"].resolve()),
           str(files["str"].resolve()), str(cfg.path(ff).resolve())]
    proc = subprocess.run(cmd, cwd=out_dir, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"cgenff_charmm2gmx failed:\n{proc.stdout}\n{proc.stderr}")
    # normalise names to <RESNAME>.itp / .prm
    for ext in ("itp", "prm"):
        low = out_dir / f"{files['resname'].lower()}.{ext}"
        if low.exists() and low != files[ext]:
            shutil.move(low, files[ext])
    return [" ".join(cmd), proc.stdout]
