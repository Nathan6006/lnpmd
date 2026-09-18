# mdp/

GROMACS run-parameter files, **generated** from `config/md.yaml` by

    python scripts/make_mdp.py            # real settings (needs all md.yaml decisions filled in)
    python scripts/make_mdp.py --smoke    # 1000-step test versions -> mdp/smoke/

The files here are for reading and review. Each run directory gets its own
copy (`runs/<system>/rep<N>/mdp/`, with that replica's velocity seed) when
`scripts/run.py` prepares it. Never edit an .mdp by hand: change md.yaml and
regenerate, so every run stays consistent.

The nonbonded block (PME, 1.2 nm, force-switch from 1.0 nm, no dispersion
correction) is hardcoded in `src/lipidmd/mdp.py` because CHARMM36 requires it.
