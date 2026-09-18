# inputs/

Files produced by the manual steps (docs/MANUAL_STEPS.md). Expected layout:

    inputs/
      charmm-gui/
        popc_pops_chol/gromacs/   step5_input.gro, topol.top, toppar/ ...
        popc/gromacs/
      cgenff/
        ECO/    ECO0.mol2 ECO0.str  ECO1.mol2 ECO1.str  ECO2.mol2 ECO2.str
        MC3/    MC3.mol2  MC3.str   MC3H.mol2 MC3H.str
        SM102/  SM12.mol2 SM12.str  SM1H.mol2 SM1H.str
      params/                    (written by scripts/prepare_params.py)
        MC3/    MC3.itp MC3.prm MC3H.itp MC3H.prm
        ...
        extra_atomtypes.itp      (only if needed, see MANUAL_STEPS troubleshooting)
      aggregates/
        ECO_aggregate.pdb  MC3_aggregate.pdb  SM102_aggregate.pdb   (64 NEUTRAL molecules each)

File names follow the residue names in config/lipids.yaml.
