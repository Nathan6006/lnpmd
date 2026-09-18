#!/usr/bin/env python3
"""A stand-in for `gmx` used only by tests/test_build.py.

Implements just enough of -version / solvate / grompp / genion for the build
pipeline to be exercised without GROMACS installed. It is NOT a simulation
engine: grompp copies coordinates into the "tpr", genion swaps waters for
ions in order.
"""
import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from lipidmd import topology as top  # noqa: E402

args = sys.argv[1:]


def opt(flag):
    return args[args.index(flag) + 1]


if args[0] == "-version":
    print(f"Data prefix:  {os.environ['FAKE_GMX_PREFIX']}")
elif args[0] == "solvate":
    solute, box = top.read_gro(opt("-cp"))
    water, wbox = top.read_gro(opt("-cs"))
    sol_xyz = np.vstack([r.xyz for r in solute])
    out = list(solute)
    reps = np.ceil(box / wbox).astype(int)
    for i in range(reps[0]):
        for j in range(reps[1]):
            for k in range(reps[2]):
                shift = wbox * [i, j, k]
                for w in water:
                    xyz = w.xyz + shift
                    if (xyz >= box).any():
                        continue
                    if np.min(np.linalg.norm(sol_xyz - xyz[0], axis=1)) < 0.3:
                        continue
                    out.append(top.Residue(w.resname, list(w.names), xyz))
    top.write_gro(opt("-o"), out, box)
elif args[0] == "grompp":
    Path(opt("-o")).write_text(Path(opt("-c")).read_text())   # "tpr" = the coordinates
elif args[0] == "genion":
    res, box = top.read_gro(opt("-s"))
    n_p, n_n = int(opt("-np")), int(opt("-nn"))
    waters = [i for i, r in enumerate(res) if r.resname == "TIP3"]
    take = waters[-(n_p + n_n):] if n_p + n_n else []
    ions = [top.Residue(opt("-pname") if k < n_p else opt("-nname"),
                        [opt("-pname") if k < n_p else opt("-nname")], res[i].xyz[:1])
            for k, i in enumerate(take)]
    keep = [r for i, r in enumerate(res) if i not in set(take)]
    top.write_gro(opt("-o"), keep + ions, box)
else:
    sys.exit(f"fake gmx: unsupported {args}")
