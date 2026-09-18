"""lipidmd: all-atom MD pipeline for ionizable lipid / membrane interaction.

Module map (see README.md for the workflow):
    config       load + sanity-check config/*.yaml
    protonation  which molecules carry which protons
    params       CGenFF .str parsing, penalty checker, GROMACS parameter prep
    topology     reading/writing GROMACS .itp/.top/.gro/.ndx files
    build        assemble a simulation box from membrane + aggregate
    mdp          generate GROMACS run-parameter (.mdp) files
    run          SLURM job scripts, checkpoint-restart chaining, status
    analysis/    one module per observable
    report       mean +/- SD table across replicas
    validation   pure-POPC check against published CHARMM36 numbers
"""

__version__ = "0.1.0"
