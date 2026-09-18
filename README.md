# lipidmd

All-atom molecular dynamics (MD) pipeline for comparing how three ionizable
lipids (**ECO**, **DLin-MC3-DMA**, **SM-102**) interact with an
endosome-like membrane as a function of **protonation state**. Lu lab, CWRU.

**Why:** ionizable lipids are near-neutral at blood pH (7.4) and become
positively charged in the acidic endosome (pH ~5–6.5). The endosomal-escape
hypothesis says the charged lipids then grab the negatively charged endosomal
membrane (POPS) and disrupt it. This pipeline simulates each lipid at 0%,
50% and 100% protonation next to a POPC/POPS/cholesterol membrane and
measures what happens.

```
 manual (web)                    this code                               cluster      Mac
 ─────────────                   ─────────                               ───────      ───
 CHARMM-GUI membrane ─┐
 CGenFF parameters  ──┼─> prepare_params ─> run.py: build box ─> job.sh ─> GROMACS ─> analyze.py ─> report.py
 aggregate .pdb     ──┘                     (protonate, solvate, ions)    (em, eq, prod)
                                            validate.py: pure POPC must match published CHARMM36 first
```

- Design choices and their reasoning: [`docs/DECISIONS.md`](docs/DECISIONS.md)
- Web-service steps, click by click: [`docs/MANUAL_STEPS.md`](docs/MANUAL_STEPS.md)

---

## A 5-minute MD primer (the words used below)

- **Force field:** the set of equations and numbers (bond stiffness,
  partial charges, van der Waals sizes) that says how atoms push and pull on
  each other. We use **CHARMM36m** for lipids/water/ions and **CGenFF** for
  our three ionizable lipids.
- **Topology (`.top`, `.itp`):** which atoms exist, how they're bonded and
  which force-field parameters apply. **Coordinates (`.gro`, `.pdb`):**
  where the atoms are.
- **`.mdp`:** the settings for one GROMACS run (timestep, temperature,
  cutoffs, ...). **`grompp`** combines mdp + coordinates + topology into a
  **`.tpr`**, and **`mdrun`** runs the `.tpr`.
- **Minimization → equilibration → production:** first remove bad contacts
  from building the box; then let water and ions relax around restrained
  lipids, releasing the restraints step by step; then run unrestrained.
  Only production is analyzed.
- **Periodic boundary conditions:** the box repeats infinitely in x, y and
  z. A molecule leaving the top re-enters at the bottom. The membrane spans
  the box in x/y, so it is effectively infinite.
- **Checkpoint (`.cpt`):** a snapshot of everything needed to continue a
  run exactly. Long runs are split into many SLURM jobs that each continue
  from the last checkpoint.
- **Replica:** an independent copy of the same system started with
  different random velocities. Differences between systems count only if
  they're bigger than the scatter between replicas.

---

## Repository layout

```
config/          all settings, one YAML file per topic (null = still to decide)
  lipids.yaml      the 3 ionizable lipids: names, formula, amine atoms, charge-state residue names
  membranes.yaml   the 2 membranes; `target:` switches between them
  systems.yaml     3x3 matrix, replicas, seeds, aggregate placement
  md.yaml          temperature, pressure, salt, thermostat/barostat, stages, lengths
  hpc.yaml         SLURM partition, GPU, wall clock, module lines (TODO)
  analysis.yaml    analysis definitions (cutoffs, selections)
  validation.yaml  pure-POPC reference values + citations
  paths.yaml       where manual-step outputs and tools live
mdp/             generated .mdp files (reference copies)
src/lipidmd/     the Python package (see src/lipidmd/__init__.py for a module map)
  analysis/        one module per observable
scripts/         command-line entry points (thin wrappers)
inputs/          manual-step outputs go here (see inputs/README.md)
validate/        validation report + PASSED marker
docs/            DECISIONS.md, MANUAL_STEPS.md
tests/           pytest suite
```

Generated at run time (not in git): `runs/`, `runs_smoke/`, `results/`.

---

## Step 1 — Install

On the Mac (and the same on the cluster, see the note at the end of this step):

```bash
git clone <this repo> lnpmd && cd lnpmd
conda env create -f environment.yml     # installs Python deps, MDAnalysis, GROMACS (CPU) and this package
conda activate lipidmd
pytest                                  # all tests should pass
```

No conda? Use `python3 -m venv .venv && source .venv/bin/activate && pip install -e ".[dev]"`
instead. For GROMACS on the Mac, `brew install gromacs` also works.

On the cluster, use the cluster's GPU build of GROMACS (`module load ...`,
recorded in `config/hpc.yaml`) instead of the conda one. Delete the
`gromacs` line from `environment.yml` there, or just
`pip install -e .` into any Python ≥ 3.10.

## Step 2 — Look at the configuration

```bash
python scripts/check_config.py
```

This checks that the config files are consistent and lists every value that
is still `null`. Nothing in this project silently picks a scientific value
for you. Each null is either:

- **TODO**: a fact you need to supply (atom names, HPC partition, ...), or
- **ASK**: a scientific decision not yet made (thermostat, equilibration
  schedule, contact cutoff, ...). Each has a comment in the YAML, and
  `docs/DECISIONS.md` section C explains what's at stake.

You can do Steps 3–4 (a smoke run) before settling the ASK items.

## Step 3 — Manual inputs

Follow [`docs/MANUAL_STEPS.md`](docs/MANUAL_STEPS.md):

1. Build both membranes on CHARMM-GUI → `inputs/charmm-gui/<membrane>/gromacs/`
2. Make one `.mol2` per lipid charge state, run CGenFF → `inputs/cgenff/<LIPID>/`
3. Fill in the amine atom names in `config/lipids.yaml`, then check and convert:
   ```bash
   python scripts/check_cgenff.py --all          # penalty report only
   python scripts/prepare_params.py --lipid MC3  # penalties + formula + conversion + name checks
   ```
4. Put a 64-molecule neutral aggregate in `inputs/aggregates/MC3_aggregate.pdb`.

To reach a smoke run fastest, do this for **one** lipid (MC3 is the
simplest: one amine) and the target membrane.

## Step 4 — Smoke run (pipeline test)

A smoke run does 1000 steps per stage. It tests that everything connects:
building, protonation, topology, grompp, mdrun, checkpointing, centring.
**It is not science.** Open ASK decisions are filled with placeholder
values, and the script prints a banner listing exactly which ones.

```bash
# 1. build the box and write every input file, submit nothing
python scripts/run.py --smoke --systems MC3_p050 --replicas 1 --dry-run
```

Look at what it made in `runs_smoke/MC3_p050/rep1/`:

| file | what it is |
|---|---|
| `start.gro` | the built box (open it in VMD or PyMOL) |
| `topol.top`, `toppar/`, `ionizable/` | topology |
| `index.ndx` | atom groups: MEMB, AGG, MEMB_AGG, SOLV |
| `build_manifest.json` | every choice the build made: which molecules are protonated, ion counts, box size |
| `mdp/*.mdp` | one settings file per stage |
| `job.sh` | the SLURM script (with `TODO_...` where hpc.yaml is unset) |
| `build.log` | GROMACS output from the build |

```bash
# 2a. run it right here on the Mac (the #SBATCH lines are just comments to bash)
bash runs_smoke/MC3_p050/rep1/job.sh

# 2b. or, on the cluster once hpc.yaml is filled in:
python scripts/run.py --smoke --systems MC3_p050 --replicas 1

# 3. watch it
python scripts/status.py --smoke
```

A successful smoke run ends with `.state` = `DONE` and a
`prod_centered.xtc`. If a stage fails, the reason is in `.failure` and the
stage's `.log`.

You can also smoke-test the analysis on it:
```bash
python scripts/analyze.py --smoke --systems MC3_p050 --replicas 1
```
Modules whose ASK settings are still null report `SKIPPED: ... (open
decision)` rather than guessing.

## Step 5 — Validation (must pass before production)

A pure POPC bilayer has to reproduce published CHARMM36 area per lipid,
thickness and order parameters. If it doesn't, something in the setup is
wrong, and nothing else can be trusted.

1. Fill in `config/validation.yaml`: temperature, length, discard time, and
   the reference values/tolerances from the cited papers. They are null on
   purpose, so copy them from the papers yourself.
2. Run:
   ```bash
   python scripts/validate.py prepare      # builds runs/validation/popc/rep1
   python scripts/validate.py submit
   python scripts/validate.py status
   # ...when DONE, copy prod.tpr + prod_centered.xtc back (see step 6), then:
   python scripts/validate.py analyze
   python scripts/validate.py compare      # prints PASS/FAIL per quantity
   ```
3. `compare` writes `validate/report.txt`, and on a full pass also
   `validate/PASSED`. **`scripts/run.py` refuses to submit production runs
   until `validate/PASSED` exists.** Delete it and re-validate whenever you
   change md.yaml, force-field files or the GROMACS version.

## Step 6 — Production runs, status, analysis, report

Once every ASK item in md.yaml/systems.yaml is decided and validation passed:

```bash
python scripts/run.py --dry-run          # build all 27, inspect
python scripts/run.py                    # build (if needed) + submit all 27
python scripts/status.py                 # per run: queued / running / done / failed / stalled + ns simulated
```

Long runs continue automatically. Each job stops cleanly shortly before the
wall-clock limit, writes a checkpoint and resubmits itself; the next job
continues exactly where it stopped (see `src/lipidmd/run.py`). If status
says **stalled** (SLURM killed the job before it could checkpoint), look at
the newest `slurm-*.out` and run `sbatch runs/<system>/rep<N>/job.sh`. It
resumes from the last checkpoint.

**Copy results to the Mac.** Only what the analysis needs. The full
trajectory stays on the cluster.
```bash
rsync -av --include='*/' --include='prod.tpr' --include='prod_centered.xtc' \
      --include='build_manifest.json' --exclude='*' \
      cluster:path/to/lnpmd/runs/ runs/
```

**Analyze and report:**
```bash
python scripts/analyze.py                 # all runs with a trajectory -> results/<system>/rep<N>/*.csv
python scripts/report.py                  # results/summary_table.{csv,md}: 9 rows, mean ± SD over replicas
```

### What each analysis writes

| module | CSV | measures |
|---|---|---|
| `membrane` | `membrane.csv` | area per lipid, bilayer thickness (P–P), per frame |
| `order` | `order.csv` | S_CD per tail carbon (POPC, POPS), per frame |
| `insertion` | `insertion.csv` | each ionizable lipid's COM and amine z relative to the membrane centre, with its charge |
| `contacts` | `contacts.csv`, `contact_events.csv` | amine N ↔ POPS phosphate/carboxylate O: counts per frame; each contact event and its duration (residence time) |
| `sasa` | `sasa.csv` | total and hydrophobic SASA of the aggregate |
| `water` | `water_density.csv`, `water_core.csv`, `water_permeation.csv` | water density profile across the membrane; waters in the core per frame; full permeation events |
| `convergence` | `convergence_running.csv`, `convergence_blocking.csv` | running averages; block-averaged standard errors |

All CSVs are "tidy": one row per observation, with `system` and `replica`
columns, ready for pandas/R/Excel.

---

## Command reference

| command | does |
|---|---|
| `scripts/check_config.py` | consistency check + list of open decisions |
| `scripts/check_cgenff.py FILE.str [--lipid L] [--threshold 10]` | CGenFF penalty report (`--all` for every configured file) |
| `scripts/prepare_params.py --lipid L \| --all` | penalties, formula/charge checks, .str → GROMACS conversion |
| `scripts/make_mdp.py [--smoke]` | write reference mdp files to `mdp/` |
| `scripts/run.py [--systems ...] [--replicas ...] [--dry-run] [--smoke] [--skip-build] [--rebuild]` | build + write inputs + submit |
| `scripts/status.py [--smoke]` | state and simulated time of every run |
| `scripts/validate.py prepare\|submit\|status\|analyze\|compare` | validation harness |
| `scripts/analyze.py [--systems ...] [--observables ...]` | analysis CSVs |
| `scripts/report.py` | the summary table |

Every script has `--help`.

## Tests

```bash
pytest            # config loading, CGenFF parser, protonation, mdp/job scripts, build (with a fake gmx), analysis
```

`tests/test_protonation.py` pins down the protonation rule (ECO 50% → all
+1; MC3 50% → half +1). `tests/test_build.py` runs the whole build against a
stub `gmx`, so it runs without GROMACS installed. The analysis tests need
MDAnalysis and are skipped without it.
