# Figures & Overall Results

Summary figures and the overall interpretation for the
[Rosavin–H₂O₂ mechanism study](../README.md), combining the MD, HOMO/LUMO
(DFT), and BDE analyses.

## Contents

```
figures/
├── FIGURES_README.md      # this file
├── make_summary_figure.py # generates results_summary.png from Table 1 & 2 values
└── results_summary.png    # MD + HOMO/LUMO summary (BDE panel pending full data)
```

## Overall results

**MD:** No significant compound-specific difference in direct H₂O₂
interaction (distance, H-bonds, SASA).

**HOMO/LUMO:** The cinnamyl group (rosavin & rosin) confers greater
electron-donating capacity than L-arabinose alone. Rosavin vs. rosin is
a small, environment-dependent effect.

**BDE:** Rosavin has more low-BDE hydroxyl sites than rosin (4/6 vs 0/4
below 110 kcal/mol) and its single most reactive site is 3.9 kcal/mol
lower than rosin's.

Taken together, the frontier-orbital and bond-dissociation-energy
analyses indicate that possession of the cinnamyl (π-conjugated) moiety,
shared by rosavin and rosin but absent in L-arabinose, confers
substantially greater intrinsic oxidizability, as reflected in a clear
and robust HOMO difference (0.75–0.85 eV in both gas phase and water)
that exceeds typical DFT uncertainty. Within this shared
cinnamyl-glycoside scaffold, rosavin possesses a greater number of
low-BDE hydroxyl sites than rosin (4 of 6 vs. 0 of 4 below 110 kcal/mol)
and its most readily accessible O–H site is modestly lower in BDE than
rosin's (by ~3.9 kcal/mol); as an exhaustive, non-selective comparison
across all hydroxyl positions in each molecule, this pattern is reported
with confidence, though it is presented as a descriptive structural
trend rather than a formally tested statistical effect.
