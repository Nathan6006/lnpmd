ECO STRUCTURE PLACEHOLDER
=========================
Put the ECO structure files here (you supply them; nothing was guessed):

  ECO0.mol2 / ECO0.str   neutral
  ECO1.mol2 / ECO1.str   +1  (proton on the site listed first in lipids.yaml protonation_order)
  ECO2.mol2 / ECO2.str   +2  (both amines protonated)

Then fill in config/lipids.yaml -> ECO: formula_neutral, mw_neutral,
protonatable_sites[*].atom / added_h, protonation_order.
See docs/MANUAL_STEPS.md step 2.
