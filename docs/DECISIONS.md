# Decisions

This file records every scientific choice in the pipeline: what was chosen,
why, and where it lives in the code or config, so each one can be defended
or changed. There are three sections:

- **A. Specified by the project brief.** Settled.
- **B. Chosen during implementation.** Technical choices that don't change
  the physics being sampled, or that follow directly from section A. Review
  them, but none were meant to be scientific.
- **C. Open questions.** Deliberately left as `null` in `config/`. Nothing
  guesses these for you. `python scripts/check_config.py` lists them.

---

## A. Specified choices

### Biological question and design
- **Question:** does protonation of an ionizable lipid (near-neutral at
  pH 7.4, cationic in the acidic endosome) let it engage and disrupt an
  anionic endosome-like membrane? This is the endosomal-escape hypothesis.
- **Design:** 3 lipids × 3 protonation levels = 9 systems, × 3 replicas = 27
  runs. Replicas are needed because a single MD trajectory only samples one
  path. Differences between systems count only if they exceed the
  replica-to-replica spread (report.py reports mean ± SD across replicas).
  → `config/systems.yaml`

### Lipids
- **ECO** (the lab's lipid): two protonatable amines, one primary and one
  tertiary. **DLin-MC3-DMA** and **SM-102**: one tertiary amine each. These
  are the benchmark lipids from approved LNP products.
- SMILES are **deliberately blank**, to be verified against PubChem.
  IUPAC names and formulas for MC3/SM-102 were counted from the IUPAC names
  and must be confirmed. Apparent pKa values (MC3 6.44, SM-102 6.68) are
  LNP TNS-assay values from the cited papers. They're informational only
  and nothing is computed from them. Note that an apparent formulation pKa
  isn't a molecular pKa. → `config/lipids.yaml`

### Protonation
- **Definition:** fraction of **all ionizable amines** that are protonated,
  so the definition is the same for one- and two-amine lipids.
  0% ≈ pH 7.4, 50% ≈ pH 6.5, 100% ≈ pH 5.0. The pH labels are nominal.
- **Distribution rule:** protons are spread over molecules as evenly as
  possible. Any remainder goes to randomly chosen molecules with a fixed
  seed. This one rule gives exactly the two required behaviours: ECO 50% →
  every molecule +1; MC3/SM-102 50% → half +1, half neutral, chosen
  randomly. It also stays well defined for any other level (e.g. ECO 75% →
  half +1, half +2).
  → `src/lipidmd/protonation.py`, tested in `tests/test_protonation.py`
- **Fixed protonation states, no constant-pH MD in v1.** This is simpler
  and fully reproducible. The limitation is that protons can't move in
  response to the local environment (e.g. a lipid buried in the membrane
  keeping its charge).
- Each charge state is a **separately parameterized molecule** (own CGenFF
  run). Protonation changes atom types and partial charges near the N, not
  just the net charge.

### System contents
- Only the raw ionizable lipid: no helper lipid, no cholesterol in the
  aggregate, no mRNA. This isolates the ionizable lipid's own behaviour.
- About 64 molecules as a pre-equilibrated aggregate placed next to the
  membrane. `n_molecules: 64` → `systems.yaml`

### Membranes
- **Target:** POPC/POPS/cholesterol 50:30:20, symmetric, 256 lipids per
  leaflet. POPS (net −1) provides the anionic surface that the
  electrostatic mechanism needs. Cholesterol gives realistic packing.
- **Control:** pure POPC, same size. Zwitterionic, so electrostatic effects
  should be much weaker.
- Switching between them is one line: `membranes.yaml → target`.

### Conditions
- 310 K, 1 bar, **semi-isotropic** pressure coupling. The membrane plane
  and the normal scale independently, so the membrane finds its own area
  per lipid instead of having it imposed.
- 150 mM NaCl **plus** counterions to neutralize. PME electrostatics require
  a neutral box, and the 150 mM keeps physiological screening.

### Force field and simulation settings
- **CHARMM36m** for lipids, **CHARMM-modified TIP3P** water, **CGenFF** for
  the three ionizable lipids. CGenFF is designed to be compatible with
  CHARMM36.
- **PME**, 1.2 nm cutoff, **Lennard-Jones force-switch from 1.0 to
  1.2 nm**. CHARMM36 lipids were parameterized with exactly this
  treatment, and using another one gives the wrong area per lipid and
  order parameters with no error message. It's hardcoded in
  `src/lipidmd/mdp.py` (`CHARMM36_NONBONDED`) and deliberately kept out of
  config so it can't be changed by accident.
- **LINCS on bonds to hydrogen, 2 fs timestep.** Constraining the fastest
  vibrations (X–H bonds) is what makes 2 fs stable.
- GROMACS on a SLURM GPU cluster. Analysis in MDAnalysis on the Mac.

### Out of scope for v1
Coarse-grained/Martini, constant-pH MD, free energies (e.g. umbrella
sampling for a water permeation PMF), mRNA, helper lipids, plotting, web
interface.

---

## B. Chosen during implementation

| Choice | Why | Where |
|---|---|---|
| **`DispCorr = no`** | Part of the CHARMM36 cutoff protocol: the lipids were fit without long-range dispersion correction. **Please confirm**: it follows from "use the CHARMM36 scheme" but wasn't stated explicitly. | `mdp.py` |
| Proton placed along the amine's lone-pair direction, N–H 0.101 nm | Only a starting geometry. Energy minimization relaxes it. | `build.py: add_proton` |
| Water removed between the two phosphate planes after solvation | `gmx solvate` fills every gap, including the bilayer's oily core, which is unphysical. Water in the headgroup region comes back during equilibration. | `build.py` step 5 |
| All CHARMM-GUI water/ions discarded; box re-solvated with CHARMM TIP3P coordinates from GROMACS `spc216.gro`, atoms renamed | The aggregate changes box height and net charge, so solvent must be redone. Any equilibrated 3-site water box is fine as *starting coordinates*. The TIP3P parameters come from CHARMM-GUI's `TIP3.itp`. | `build.py` |
| Validation box keeps CHARMM-GUI's water thickness | No aggregate, so there's no reason to change it. | `build.py` |
| Aggregate centred over the box in x/y | Symmetric start, as far as possible from its own periodic images. | `build.py: place_aggregate` |
| EM: steepest descent, `emtol` 1000 kJ/mol/nm, max 50 000 steps | Standard way to remove clashes before dynamics. Doesn't affect sampling. | `mdp.py` |
| `nstlist` 20 with automatic Verlet buffer, `nstcalcenergy` 100 | Performance settings. GROMACS sizes the buffer to keep energy drift below its tolerance. | `mdp.py` |
| Restraint reference = the built structure for every stage (`-r start.gro`), `refcoord-scaling = com` | Restraints hold atoms at their built positions instead of drifting stage to stage. `com` is required when restraints are combined with pressure coupling. | job.sh / `mdp.py` |
| Each equilibration stage continues from the previous stage's checkpoint (`grompp -t`) | Carries exact velocities forward. | job.sh |
| Heavy-atom position restraints on the aggregate, switched by `POSRES_AGG` | Lets the aggregate be held during early equilibration. The strengths are an open decision (C). | `topology.write_posre_itp` |
| Checkpoint-restart via `mdrun -maxh` + self-resubmission; `max_segments` 50 | Long runs continue across wall-clock limits. The cap prevents endless resubmission if something is wrong. | `run.py` |
| Trajectory post-processing `trjconv -pbc mol -center` on the membrane | Analysis needs whole molecules and a centred membrane. | job.sh |
| Membrane centre = centre of mass of all membrane lipids; z-distances wrapped to [−Lz/2, Lz/2) | Standard. The wrapping handles an aggregate that drifts through the periodic boundary toward the other leaflet. | `analysis/common.py` |
| Leaflet assignment per frame by headgroup atom (P, or O3 for cholesterol) above/below the centre | Simple, and handles cholesterol flip-flop. | `analysis/membrane.py` |
| **Thickness = D_PP** (mean P–P distance between leaflets) | The most common simulation definition. The validation references must use the same definition. | `analysis.yaml` |
| S_CD from all C–H bonds (found via bonds in the .tpr), membrane normal = z | Standard NMR-comparable definition. Independent of hydrogen naming. | `analysis/order.py` |
| POPS "phosphate and carboxylate oxygens" = O11–O14 + O13A/O13B | Literal reading: all four O on P, plus both serine carboxylate O. **Verify the names in `toppar/POPS.itp`**, and decide whether the ester oxygens O11/O12 should count. | `analysis.yaml` |
| Contact = any listed POPS O within cutoff of an amine N; per (N, POPS molecule) pair | Counts one lipid–lipid contact once, however many of its O are close. | `analysis/contacts.py` |
| Residence time = event duration; events still open at the end are flagged `censored` | Censored events make the mean residence time an underestimate. That's reported, not hidden. | `analysis/contacts.py`, `report.py` |
| SASA by Shrake–Rupley, Bondi radii, 240 points per atom | Standard algorithm, implemented here to avoid another dependency. Accuracy ~1% at 240 points. | `analysis/sasa.py` |
| Water density: 0.1 nm bins, number density of water O | Resolution only. | `analysis.yaml` |
| Permeation = bulk → core → opposite bulk, without returning | Standard event definition. PBC wrap-around isn't counted. The region boundaries are open (C). | `analysis/water.py` |
| Convergence: running average + Flyvbjerg–Petersen block averaging (block sizes 1, 2, 4, …) | Parameter-free. The plateau of SEM vs block size is the error bar. | `analysis/convergence.py` |
| Report: per-replica time average after `discard_ns`, then mean ± SD over replicas | Replica SD is the honest uncertainty for comparing systems. | `report.py` |
| Replica seeds: fixed list `[1101, 2202, 3303]`, protonation seed 20260918, ion seed 4404 | Arbitrary, but fixed for reproducibility. | `systems.yaml` |
| Residue names ECO0/ECO1/ECO2, MC3/MC3H, SM12/SM1H | Up to 4 characters, so they fit PDB format. Change them freely, but keep them consistent with the file names. | `lipids.yaml` |
| **Validation must pass before production submits** | Brief: "must pass before I trust anything else". Enforced with `validate/PASSED`. | `run.py`, `validation.py` |

### CGenFF penalties
Policy built into `check_cgenff.py`, using the CGenFF authors' thresholds:
report everything > 10 and flag > 50 as needing manual attention. Anything
within 3 bonds of an ionizable N is marked, because those parameters
control exactly the chemistry being compared across protonation states. How
to handle high penalties near the N (e.g. refit dihedrals with the FFParam
tool, or compare against a validated analog) is a decision for you and
Dr. Lu once the real `.str` files exist.

---

## C. Open questions

Each item below is `null` in `config/` until you decide it. Where there is a
common choice I've mentioned it as a starting point for the discussion, not
as a recommendation. `--smoke` runs fill these with placeholders
(`SMOKE_FALLBACKS` in `src/lipidmd/config.py`), print a banner saying so,
and write to `runs_smoke/`. Those placeholders are **not** decisions.

### Membrane
1. **Exact per-leaflet counts** for 50:30:20 of 256: 76.8 POPS and 51.2
   CHL1 aren't integers. One option is 128/77/51. → `membranes.yaml`
2. **APL with cholesterol:** does cholesterol count in "lipids per leaflet"?
   Conventions differ, so the choice decides which literature values are
   comparable. → `analysis.yaml`

### System assembly
3. **Aggregate preparation:** how are the 64-molecule aggregates made and
   pre-equilibrated? Is one neutral aggregate used for all protonation levels
   (protons added afterwards, as the build currently does), or is each level
   pre-equilibrated separately? → `MANUAL_STEPS.md` step 4
4. **Initial gap** between the aggregate and the membrane surface, and
   **water padding**. A gap too small biases toward contact, and one too
   large wastes simulation time before anything happens. The padding must
   exceed the 1.2 nm cutoff with margin, so the aggregate doesn't interact
   with the membrane's periodic image. → `systems.yaml → aggregate`
5. **ECO +1 site:** at 50%, is the proton on the primary or the tertiary
   amine? ECO's microscopic pKa values (measured or computed) are the
   natural way to settle it. → `lipids.yaml → ECO.protonation_order`
6. **Replica variation:** velocities only, or also a different random choice
   of which MC3/SM-102 molecules are charged? → `systems.yaml`
7. **Salt basis:** is 150 mM relative to the water volume (what the water
   actually sees) or to the whole box (`gmx genion -conc`, which counts
   the membrane volume and adds more ions)? → `md.yaml → salt_basis`

### MD protocol
8. **Thermostat** and τ_T. V-rescale, τ_T = 1 ps is what CHARMM-GUI
   generates.
9. **Coupling groups** for the thermostat and centre-of-mass removal.
   CHARMM-GUI uses membrane / solvent. The question is where the aggregate
   goes: with the membrane (`MEMB_AGG`, `SOLV`), with the solvent, or as its
   own group. Removing COM motion per group affects how the aggregate may
   drift relative to the membrane.
10. **Barostat** (equilibration and production), τ_P, compressibility.
    C-rescale is correct for both in modern GROMACS. Parrinello–Rahman is
    also common for production. The compressibility of water is 4.5e-5 bar⁻¹.
11. **Equilibration schedule:** number of stages, NVT/NPT, lengths, and
    restraint strengths on lipid headgroups, lipid dihedrals and the
    aggregate. CHARMM-GUI's `step6.*.mdp` files are a published reference
    schedule. Note that its first stages use 1 fs, while the brief fixes
    2 fs.
12. **Minimization restraints.**
13. **Production length.** The brief mentions 500 ns as an example. The
    convergence output from a first run should inform this.
14. **Output intervals** (xtc, energy, log). These set the time resolution
    of residence times and permeation detection, and the disk use.

### Analysis
15. **Equilibration discard** (`discard_ns`). Decide after looking at the
    running averages from a first run.
16. **Contact cutoff** and **gap tolerance** for residence times.
17. **SASA:** probe radius (0.14 nm is the usual water probe), which atoms
    count as hydrophobic, and whether membrane atoms bury aggregate
    surface (SASA "in context" vs the aggregate's own surface).
18. **Permeation region boundaries:** core half-width and bulk boundary.
    Looking at the water density profile from a first run is the natural way
    to set these.

### Validation
19. **Temperature** (published CHARMM36 POPC data are mostly at 303 K),
    length, discard time.
20. **Reference values and tolerances.** Each must be taken from the cited
    paper's tables with the same definitions as ours (D_PP vs D_HH; signed
    S_CD). None were filled in from memory.

### Infrastructure (facts, not decisions)
21. HPC partition, GPU type, CPUs per GPU, wall-clock limit, module lines.
    → `hpc.yaml`
22. Paths to `cgenff_charmm2gmx` and the CHARMM36 GROMACS port.
    → `paths.yaml`
23. Atom names of the protonatable N and the added H for each lipid. These
    come from your structure files. → `lipids.yaml`
