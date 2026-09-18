"""Small readers/writers for GROMACS text formats.

Kept dependency-free (no MDAnalysis) so the build step runs anywhere GROMACS
runs. Only the parts of each format this project needs are handled.

Formats in one sentence each:
    .gro  coordinates in nm, fixed-width columns, box vector on the last line
    .pdb  coordinates in Angstrom (we convert to nm on read)
    .itp  "include topology": one or more [ moleculetype ] definitions
    .top  the full topology: includes + [ system ] + [ molecules ] counts
    .ndx  named groups of 1-based atom numbers ("index groups")
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np


# --------------------------------------------------------------------------
# Coordinates
# --------------------------------------------------------------------------
@dataclass
class Residue:
    """One molecule's worth of atoms (every molecule here is one residue)."""

    resname: str
    names: list[str]
    xyz: np.ndarray  # (n_atoms, 3), nm

    def copy(self) -> "Residue":
        return Residue(self.resname, list(self.names), self.xyz.copy())

    def index(self, atom_name: str) -> int:
        try:
            return self.names.index(atom_name)
        except ValueError:
            raise KeyError(f"atom {atom_name!r} not in residue {self.resname}") from None


def read_gro(path: str | Path) -> tuple[list[Residue], np.ndarray]:
    """Residues and box (3,) in nm. Only rectangular boxes are supported."""
    lines = Path(path).read_text().splitlines()
    n = int(lines[1].strip())
    residues: list[Residue] = []
    key = None
    names: list[str] = []
    coords: list[list[float]] = []
    for line in lines[2 : 2 + n]:
        this_key = (line[0:5], line[5:10].strip())
        if this_key != key and names:
            residues.append(Residue(key[1], names, np.array(coords)))
            names, coords = [], []
        key = this_key
        names.append(line[10:15].strip())
        coords.append([float(line[20:28]), float(line[28:36]), float(line[36:44])])
    if names:
        residues.append(Residue(key[1], names, np.array(coords)))
    box_fields = [float(x) for x in lines[2 + n].split()]
    if len(box_fields) > 3 and any(abs(x) > 1e-6 for x in box_fields[3:]):
        raise ValueError(f"{path}: triclinic boxes are not supported")
    return residues, np.array(box_fields[:3])


def write_gro(path: str | Path, residues: list[Residue], box: np.ndarray, title: str = "lipidmd") -> None:
    out = [title, f"{sum(len(r.names) for r in residues):5d}"]
    atom_nr = 0
    for resnr, res in enumerate(residues, start=1):
        for name, (x, y, z) in zip(res.names, res.xyz):
            atom_nr += 1
            # .gro columns wrap at 5 digits; GROMACS itself does the same.
            out.append(f"{resnr % 100000:5d}{res.resname:<5.5s}{name:>5.5s}{atom_nr % 100000:5d}"
                       f"{x:8.3f}{y:8.3f}{z:8.3f}")
    out.append(f"{box[0]:10.5f}{box[1]:10.5f}{box[2]:10.5f}")
    Path(path).write_text("\n".join(out) + "\n")


def read_pdb(path: str | Path) -> list[Residue]:
    """ATOM/HETATM records -> residues (nm). A new residue starts whenever
    chain, residue number or residue name changes."""
    residues: list[Residue] = []
    key = None
    names: list[str] = []
    coords: list[list[float]] = []
    for line in Path(path).read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        this_key = (line[21], line[22:27], line[17:21].strip())
        if this_key != key and names:
            residues.append(Residue(key[2], names, np.array(coords)))
            names, coords = [], []
        key = this_key
        names.append(line[12:16].strip())
        coords.append([float(line[30:38]) / 10, float(line[38:46]) / 10, float(line[46:54]) / 10])
    if names:
        residues.append(Residue(key[2], names, np.array(coords)))
    return residues


# --------------------------------------------------------------------------
# Topology (.itp / .top)
# --------------------------------------------------------------------------
@dataclass
class MoleculeType:
    name: str
    atom_names: list[str] = field(default_factory=list)
    atom_types: list[str] = field(default_factory=list)
    charges: list[float] = field(default_factory=list)
    masses: list[float] = field(default_factory=list)
    resnames: list[str] = field(default_factory=list)
    bonds: list[tuple[int, int]] = field(default_factory=list)  # 0-based atom indices

    @property
    def charge(self) -> float:
        return float(sum(self.charges))

    def neighbors(self, atom_name: str) -> list[str]:
        i = self.atom_names.index(atom_name)
        out = []
        for a, b in self.bonds:
            if a == i:
                out.append(self.atom_names[b])
            elif b == i:
                out.append(self.atom_names[a])
        return out


_SECTION = re.compile(r"^\s*\[\s*([a-z_]+)\s*\]")


def _clean(line: str) -> str:
    return line.split(";", 1)[0].strip()


def parse_itp(path: str | Path) -> dict[str, MoleculeType]:
    """All moleculetypes defined in one .itp (includes are NOT followed).

    Lines inside #ifdef blocks are read like any other (they only contain
    restraints in the files we use, which we don't parse)."""
    mols: dict[str, MoleculeType] = {}
    section = None
    current: MoleculeType | None = None
    for raw in Path(path).read_text().splitlines():
        m = _SECTION.match(raw)
        if m:
            section = m.group(1)
            continue
        line = _clean(raw)
        if not line or line.startswith("#"):
            continue
        f = line.split()
        if section == "moleculetype":
            current = MoleculeType(f[0])
            mols[current.name] = current
        elif section == "atoms" and current is not None:
            # nr type resnr residue atom cgnr charge [mass]
            current.atom_types.append(f[1])
            current.resnames.append(f[3])
            current.atom_names.append(f[4])
            current.charges.append(float(f[6]))
            current.masses.append(float(f[7]) if len(f) > 7 else float("nan"))
        elif section == "bonds" and current is not None:
            current.bonds.append((int(f[0]) - 1, int(f[1]) - 1))
    return mols


def parse_atomtype_names(path: str | Path) -> set[str]:
    """Names of every [ atomtypes ] entry in a file."""
    names = set()
    section = None
    for raw in Path(path).read_text().splitlines():
        m = _SECTION.match(raw)
        if m:
            section = m.group(1)
            continue
        line = _clean(raw)
        if section == "atomtypes" and line and not line.startswith("#"):
            names.add(line.split()[0])
    return names


def parse_top_includes(path: str | Path) -> list[str]:
    """The #include "..." targets of a .top, in order."""
    out = []
    for raw in Path(path).read_text().splitlines():
        m = re.match(r'^\s*#include\s+"([^"]+)"', raw)
        if m:
            out.append(m.group(1))
    return out


def run_length(residues: list[Residue], resname_to_moltype: dict[str, str]) -> list[tuple[str, int]]:
    """[ molecules ] entries: consecutive identical molecules collapsed.
    GROMACS requires this list to match the coordinate order exactly."""
    out: list[tuple[str, int]] = []
    for res in residues:
        try:
            mt = resname_to_moltype[res.resname]
        except KeyError:
            raise KeyError(f"no moleculetype for residue {res.resname!r}") from None
        if out and out[-1][0] == mt:
            out[-1] = (mt, out[-1][1] + 1)
        else:
            out.append((mt, 1))
    return out


def write_top(path: str | Path, includes: list[str], molecules: list[tuple[str, int]],
              title: str = "lipidmd system") -> None:
    lines = ["; Generated by lipidmd.build -- regenerate rather than editing by hand.", ""]
    lines += [f'#include "{inc}"' for inc in includes]
    lines += ["", "[ system ]", title, "", "[ molecules ]", "; name      count"]
    lines += [f"{name:<10s} {count:6d}" for name, count in molecules]
    Path(path).write_text("\n".join(lines) + "\n")


def write_ndx(path: str | Path, groups: dict[str, list[int]]) -> None:
    """groups: name -> 0-based atom indices (written 1-based, as GROMACS wants)."""
    lines = []
    for name, idx in groups.items():
        lines.append(f"[ {name} ]")
        nums = [str(i + 1) for i in idx]
        for k in range(0, len(nums), 15):
            lines.append(" ".join(nums[k : k + 15]))
        lines.append("")
    Path(path).write_text("\n".join(lines))


def write_posre_itp(path: str | Path, mol: MoleculeType, switch: str, fc_macro: str) -> None:
    """Position restraints on heavy atoms, active only when `switch` is defined.

    The force constant is a preprocessor macro, so each mdp stage can set its
    own strength with define = -Dswitch -Dfc_macro=VALUE."""
    lines = [f"; heavy-atom position restraints for {mol.name}", f"#ifdef {switch}",
             "[ position_restraints ]", ";  ai  funct  fcx  fcy  fcz"]
    for i, (name, mass) in enumerate(zip(mol.atom_names, mol.masses)):
        is_h = mass < 1.5 if not np.isnan(mass) else name.upper().startswith("H")
        if not is_h:
            lines.append(f"{i + 1:6d}  1  {fc_macro} {fc_macro} {fc_macro}")
    lines += ["#endif", ""]
    Path(path).write_text("\n".join(lines))


# --------------------------------------------------------------------------
# Merging GROMACS parameter (.prm) files from several CGenFF conversions
# --------------------------------------------------------------------------
_N_TYPES = {"atomtypes": 1, "bondtypes": 2, "constrainttypes": 2, "pairtypes": 2,
            "angletypes": 3, "dihedraltypes": 4, "nonbond_params": 2}


class ParameterConflict(ValueError):
    pass


def merge_prm(paths: list[str | Path]) -> str:
    """Combine several .prm files into one, dropping exact duplicates.

    Why this exists: the neutral and protonated forms of a lipid are
    parameterized separately, and both .prm files define some of the same
    bonded types. GROMACS rejects or warns on repeated definitions, so we keep
    one copy -- and stop with an error if two files give DIFFERENT values for
    the same atom types, because then one of them would be silently wrong.
    """
    # section -> key -> tuple of parameter rows (multi-row for dihedral funct 9)
    merged: dict[str, dict[tuple, tuple]] = {}
    origin: dict[tuple, str] = {}
    order: list[tuple[str, tuple]] = []
    for path in paths:
        local: dict[tuple[str, tuple], list[tuple]] = {}
        section = None
        for raw in Path(path).read_text().splitlines():
            m = _SECTION.match(raw)
            if m:
                section = m.group(1)
                continue
            line = _clean(raw)
            if not line or line.startswith("#") or section not in _N_TYPES:
                continue
            f = line.split()
            n = _N_TYPES[section]
            types = tuple(f[:n])
            if n > 1 and types[::-1] < types:
                types = types[::-1]  # a-b-c and c-b-a are the same parameter
            funct = f[n] if section != "atomtypes" else ""
            key = (types, funct)
            params = tuple(_num(x) for x in f[n + (0 if section == "atomtypes" else 1):])
            local.setdefault((section, key), []).append(params)
        for (section, key), rows in local.items():
            rows_t = tuple(rows)
            sect = merged.setdefault(section, {})
            if key in sect:
                if sect[key] != rows_t:
                    raise ParameterConflict(
                        f"[ {section} ] {' '.join(key[0])} differs between "
                        f"{origin[(section,) + key]} and {path}: {sect[key]} vs {rows_t}")
                continue
            sect[key] = rows_t
            origin[(section,) + key] = str(path)
            order.append((section, key))

    out = ["; merged by lipidmd.topology.merge_prm"]
    for section in _N_TYPES:
        keys = [k for s, k in order if s == section]
        if not keys:
            continue
        out.append(f"\n[ {section} ]")
        for key in keys:
            types, funct = key
            for row in merged[section][key]:
                out.append("  ".join([*types, *([funct] if funct else []), *(_fmt(x) for x in row)]))
    return "\n".join(out) + "\n"


def _num(x: str):
    try:
        return round(float(x), 8)
    except ValueError:
        return x


def _fmt(x) -> str:
    return f"{x:.8g}" if isinstance(x, float) else str(x)
