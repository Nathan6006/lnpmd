# Manual steps (v1)

Two inputs come from web services that we deliberately do not script:
**CHARMM-GUI** (builds the membrane) and the **CGenFF server** (parameters for
the ionizable lipids). A third input, the **pre-equilibrated aggregate**, is
still an open decision (step 4).

Website layouts change. If a button below is named slightly differently, pick
the closest match and **write down what you chose** in `docs/DECISIONS.md`.

---

## Step 1 — Membranes from CHARMM-GUI Membrane Builder

You need two membranes: `popc_pops_chol` (the target) and `popc` (control +
validation). Repeat these steps for each.

Before you start, settle the exact per-leaflet counts for the mixed membrane:
256 × 30% = 76.8 POPS isn't a whole number, so `membranes.yaml` →
`counts_per_leaflet` has to be chosen (one option: POPC 128 / POPS 77 /
CHL1 51). Enter the same numbers on the website and in `membranes.yaml`; the
build step checks that they agree.

1. Go to <https://www.charmm-gui.org> and log in. A free academic account is
   enough.
2. **Input Generator → Membrane Builder → Bilayer Builder**.
3. Choose **Membrane Only System**, since there's no protein.
4. **Box / system size:**
   - Box type: **Rectangular**. Our code only handles rectangular boxes.
   - Length of XY: choose **Numbers of lipid components**. Our lipid
     counts are fixed, and this makes the box size follow from them.
   - Water thickness: leave the default and **record the value**. It only
     matters for the validation (pure POPC) box. Production boxes are
     re-solvated by our build step.
5. **Lipid composition**, same numbers in the upper and lower leaflet
   (the membrane is symmetric):
   - `popc_pops_chol`: PC lipids → **POPC**; PS lipids → **POPS**;
     Sterols → **CHL1** (cholesterol). Use the counts from `membranes.yaml`.
   - `popc`: **POPC** only, 256 per leaflet.
   - Click **Show the system info**. Check that the net charge is −(2 × POPS
     count), then click **Next**.
6. **Ions:** check **Include Ions**, type **NaCl**, 0.15 M.
   Our build removes CHARMM-GUI's water and ions and adds its own, but this
   option makes CHARMM-GUI include the `SOD`/`CLA` topology files we need.
7. Keep clicking **Next** through the assembly steps. Each takes a minute or
   two. Save the job ID so you can come back to it.
8. **Input generation options:**
   - Force field: **CHARMM36m**
   - Input: check **GROMACS**
   - Temperature: 310 K for `popc_pops_chol`. For `popc`, use
     `validation.yaml → temperature_K` once you've decided it. It only
     affects CHARMM-GUI's own example mdp files, which we don't use.
9. Download the `.tgz`. Unpack it:
   ```
   tar xzf charmm-gui.tgz
   ```
   Then copy its **`gromacs/`** folder here:
   ```
   inputs/charmm-gui/popc_pops_chol/gromacs/     (and inputs/charmm-gui/popc/gromacs/)
   ```
   Check that it contains `step5_input.gro`, `topol.top` and `toppar/`.
10. Check the restraint macro names the build relies on:
    ```
    grep -h "ifdef\|FC" inputs/charmm-gui/popc/gromacs/toppar/POPC.itp | sort -u
    ```
    You should see `POSRES`, `POSRES_FC_LIPID`, `DIHRES`, `DIHRES_FC`. If
    they differ, copy the names into `md.yaml → restraint_macros`.
11. Worth reading, not used directly: `gromacs/step6.*_equilibration.mdp` is
    CHARMM-GUI's recommended equilibration schedule. It's a good reference
    for the `md.yaml → equilibration` decision.

## Step 2 — CGenFF parameters for the ionizable lipids

Every **charge state** of every lipid is a separate molecule for CGenFF:

| lipid | charge 0 | +1 | +2 |
|---|---|---|---|
| ECO | `ECO0` | `ECO1` | `ECO2` |
| MC3 | `MC3` | `MC3H` | — |
| SM-102 | `SM12` | `SM1H` | — |

These names are the residue names in `config/lipids.yaml`. They must be
used exactly as written, because the code finds files by these names.

### 2a. Make the neutral structure (.mol2)

1. Get a 3D structure with **all hydrogens explicit**. For MC3 and SM-102
   you can build it from the SMILES once you have verified the SMILES
   yourself. For **ECO, use the structure file you supply.**
   Avogadro works for this (**Build → Add Hydrogens**, then
   **Extensions → Optimize Geometry**).
2. Give **every atom a unique name** (N1, C1, C2, …, H1, H2, …). CGenFF and
   GROMACS both identify atoms by name.
3. Set the molecule name (the line after `@<TRIPOS>MOLECULE`) and the
   residue name in the atom lines to the charge-0 name (e.g. `MC3`).
4. Save the file as `inputs/cgenff/<LIPID>/<resname>.mol2`, e.g.
   `inputs/cgenff/MC3/MC3.mol2`.
5. In `config/lipids.yaml`, fill in `protonatable_sites[*].atom` with the
   amine nitrogen name(s) from this file.

### 2b. Make each protonated structure

1. Open the neutral `.mol2` and **add one hydrogen to the amine nitrogen**
   (for ECO +1: the site listed first in `protonation_order`; for ECO +2:
   both).
2. Name the new hydrogen as you set it in `lipids.yaml → added_h`
   (e.g. `HN1`). **Do not rename any other atom.** The build adds this proton
   to the neutral aggregate by name, and `prepare_params.py` checks that the
   names line up.
3. Set the molecule and residue name to the charged name (e.g. `MC3H`).
   Save it as `inputs/cgenff/MC3/MC3H.mol2`.

### 2c. Run CGenFF

1. Go to the CGenFF server (now hosted by SilcsBio: <https://cgenff.silcsbio.com>)
   and log in. A free academic account is enough.
2. Upload one `.mol2` at a time. Keep the default options.
3. Download the resulting stream file and save it next to the `.mol2` with
   the same base name: `inputs/cgenff/MC3/MC3H.str`.
4. Check the penalties straight away:
   ```
   python scripts/check_cgenff.py inputs/cgenff/MC3/MC3H.str --lipid MC3
   ```
   Anything marked `*` and `ATTENTION` touches the ionizable headgroup, which
   is the chemistry this project measures. Discuss it before you go further
   (see `DECISIONS.md`, "CGenFF penalties").

## Step 3 — Convert CGenFF → GROMACS

GROMACS can't read `.str` files. The MacKerell lab (who maintain CHARMM)
provide a converter:

1. From <https://mackerell.umaryland.edu/charmm_ff.shtml#gromacs> download:
   - the **CHARMM36 force field for GROMACS** (a folder like
     `charmm36-jul2022.ff`)
   - the conversion script **`cgenff_charmm2gmx_py3_nx2.py`**. It needs
     Python 3 and `networkx`, which is already in `environment.yml`.
2. Put both under `tools/` (e.g. `tools/charmm36-jul2022.ff/`) and set
   their paths in `config/paths.yaml`.
3. Run:
   ```
   python scripts/prepare_params.py --lipid MC3      # or --all
   ```
   This runs the penalty report, checks the formula and net charge, converts
   every charge state into `inputs/params/<LIPID>/<resname>.itp/.prm`, and
   checks that each protonated `.itp` equals the neutral one plus the added
   hydrogens.

## Step 4 — Pre-equilibrated aggregate (OPEN DECISION)

The build expects `inputs/aggregates/<LIPID>_aggregate.pdb`: 64 **neutral**
molecules, whose atom names match the neutral `.mol2`. The build step adds the
protons itself.

**How** to make and pre-equilibrate this aggregate hasn't been decided yet
(see `DECISIONS.md`, open questions). Some possibilities: packing 64
molecules with Packmol and running a short MD in water, or starting from a
published aggregate/LNP-core structure. There's also the question of whether
one neutral aggregate should serve all three protonation levels. Until this is
settled, use any reasonable 64-molecule cluster for **smoke tests only**.

## Step 5 — Where things run

- **Building** (`scripts/run.py`, `scripts/validate.py prepare`) calls `gmx`.
  Run it on the cluster login node after `module load` of GROMACS, or on
  the Mac with the conda `gromacs` package. The resulting run directories are
  self-contained.
- **Simulation** runs on the cluster through SLURM.
- **Analysis** runs on the Mac. Copy back only what it needs (README step 6).

---

## Troubleshooting

**"uses atom types [...] that the CHARMM-GUI force field files do not define"**
CHARMM-GUI's `toppar/forcefield.itp` may not define every CGenFF atom type
your lipid uses. Open `tools/charmm36-jul2022.ff/ffnonbonded.itp` and copy
the lines for the missing types from its `[ atomtypes ]` section into a new
file `inputs/params/extra_atomtypes.itp`, under a `[ atomtypes ]` header. If
the same file has `[ pairtypes ]` lines involving those types, copy them too
under a `[ pairtypes ]` header. The build includes this file automatically.
First check that the copied lines don't duplicate anything already in
`forcefield.itp`.

**"restraint macros [...] do not appear in the CHARMM-GUI lipid topologies"**
See step 1.10.

**`ParameterConflict` while merging .prm files**
The neutral and protonated parameterizations gave *different* values for the
same bonded type. Don't pick one arbitrarily. Look at the two `.str` files
and decide with your PI which is right.
