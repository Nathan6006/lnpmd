# validate/

Output of the validation harness (`python scripts/validate.py ...`, README step 5):

- `report.txt` — latest comparison of the pure-POPC run against `config/validation.yaml`
- `PASSED` — exists only if the latest comparison passed. `scripts/run.py`
  will not submit production runs without it.

Delete `PASSED` whenever you change anything that could affect the physics
(md.yaml, force field files, GROMACS version) and re-validate.
