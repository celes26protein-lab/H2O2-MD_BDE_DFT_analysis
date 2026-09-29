# BDE: O–H Bond Dissociation Energy

Computes the O–H bond dissociation energy (BDE) at **every** hydroxyl
site in rosavin and rosin (an exhaustive, unbiased comparison, not a
hand-picked subset), to test whether rosavin has more easily-abstracted
hydroxyl hydrogens than rosin, as part of the
[Rosavin–H₂O₂ mechanism study](../README.md).

L-arabinose is not included here — BDE was computed for rosavin and
rosin only.

## Method

For every hydroxyl (O–H) site in a compound:

1. **Parent molecule:** built from SMILES with [RDKit](https://www.rdkit.org/)
   (`EmbedMolecule` + `MMFFOptimizeMolecule`). Single-point DFT energy
   (RHF reference, closed shell) — no DFT-level re-optimization.
2. **Radical:** that one hydroxyl H atom is removed (`RWMol`), a radical
   electron is placed on its O (`SetNumRadicalElectrons(1)`,
   `SetNoImplicit(True)`), and the resulting doublet radical is given its
   own fresh 3D conformer (re-embedded and MMFF-optimized — this is a
   different structure from the parent, not just "parent minus one atom
   frozen in place"). Single-point DFT energy (UHF reference).
3. **H atom reference:** a lone H atom, doublet, UHF, single-point DFT
   energy — computed separately for gas phase and water (PCM), since
   solvation shifts it too.
4. **BDE:**

   ```
   BDE = E(radical) + E(H atom) − E(parent)
   ```

All energies: [Psi4](https://psicode.org/) (v1.11), B3LYP/6-31G*.
Implicit water uses the same IEFPCM/GePol settings as the
[HOMO/LUMO calculation](../dft_homo_lumo/DFT_README.md) (UFF radii,
`Scaling = False`, `Area = 0.3`).

## Contents

```
bde/
├── BDE_README.md   # this file
├── run_bde.py       # parent -> per-site radicals -> BDE, for one compound/environment
└── *.log            # Psi4 output logs (H atom, each parent, each radical)
```

**Not included:** `bde_results.csv` (and the equivalent HOMO/LUMO results
CSV) — omitted due to file size. The completed BDE values are reported
below and in the main README instead.

## Usage

```bash
python3 run_bde.py <name> "<SMILES>" <gas|water>
```

One run computes the H-atom reference, the parent, and **all** hydroxyl-
site radicals for that compound/environment automatically (6 sites for
rosavin, 4 for rosin):

```bash
python3 run_bde.py rosavin "O[C@H]1CO[C@@H](OC[C@H]2O[C@@H](OC\C=C\C3=CC=CC=C3)[C@H](O)[C@@H](O)[C@@H]2O)[C@H](O)[C@H]1O" gas
python3 run_bde.py rosavin "O[C@H]1CO[C@@H](OC[C@H]2O[C@@H](OC\C=C\C3=CC=CC=C3)[C@H](O)[C@@H](O)[C@@H]2O)[C@H](O)[C@H]1O" water

python3 run_bde.py rosin "C1=CC=C(C=C1)/C=C/CO[C@H]2[C@@H]([C@H]([C@@H]([C@H](O2)CO)O)O)O" gas
python3 run_bde.py rosin "C1=CC=C(C=C1)/C=C/CO[C@H]2[C@@H]([C@H]([C@@H]([C@H](O2)CO)O)O)O" water
```

## Requirements

- [Psi4](https://psicode.org/) (conda-only: `conda install -c conda-forge -c psi4 psi4`)
- [RDKit](https://www.rdkit.org/)

## Results (confirmed values)

| Compound | Site | Environment | BDE (kcal/mol) |
|---|---|---|---|
| Rosavin | O0–H30 | gas | 108.07 |
| Rosavin | O21(O23 in solvated run)–H | water | 116.69 |

Only these two site/environment combinations have been directly
confirmed against notebook output so far. Fill in the rest of the 6
(rosavin) + 4 (rosin) sites × 2 environments = 20 total values here once
`run_bde.py` has been run for all of them (or once the remaining
notebook cells are located).

Per the main README, the overall finding is: within the shared
cinnamyl-glycoside scaffold, rosavin has more low-BDE hydroxyl sites than
rosin (4 of 6 vs. 0 of 4 below 110 kcal/mol), and rosavin's most reactive
site is ~3.9 kcal/mol lower than rosin's most reactive site.

## Caveats

- **No DFT-level geometry optimization** for the parent (same limitation
  as the HOMO/LUMO calculation).
- **The radical *is* re-embedded and MMFF-optimized** as a fresh
  structure (unlike the parent, which reuses the same MMFF conformer
  throughout) — this is a meaningful methodological difference worth
  stating explicitly in the paper's Methods section if not already
  there.
- Single conformer per site (no conformational sampling of the radical).
- One site (`H30`→`O29`, gas phase) showed the farthest O···H distance in
  the intramolecular-distance screen — worth double-checking the units
  used in that screen (the printed values, labeled "nm", look more
  consistent with Å; if so, the discrepancy is only a print-label bug and
  does not affect the BDE values above).
