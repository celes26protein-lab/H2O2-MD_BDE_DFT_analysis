# Rosavin–H₂O₂ Mechanism

Computational analysis of the antioxidant mechanism of rosavin (a *Rhodiola
rosea* constituent), as described in our published paper (see
[Citation](#citation)).

## Background

*Rhodiola rosea* extract reduces hydrogen peroxide (H₂O₂) levels, and among
its major constituents, rosavin significantly decreases H₂O₂ while
salidroside has no effect. The intact phenylpropanoid glycoside structure
of rosavin is required for activity — its structural components, rosin
and L-arabinose, show no inhibitory effect on their own.

This repository investigates *why*, computationally, using a three-stage
approach:

1. **Molecular dynamics** — does rosavin physically associate with H₂O₂
   more than its structural relatives?
2. **Frontier molecular orbitals (DFT)** — is rosavin intrinsically more
   prone to donating an electron (more oxidizable)?
3. **O–H bond dissociation energy (BDE)** — how easily does each hydroxyl
   site release a hydrogen atom?

## Contents

```
md/                  # GROMACS MD setup and summary results
├── mdp/             # em.mdp, nvt.mdp, npt.mdp, md.mdp
└── topology/         # topology files, per compound
dft_homo_lumo/        # Psi4: HOMO/LUMO energies (gas phase & implicit water)
├── run_homo_lumo.py  # SMILES -> RDKit 3D geometry -> Psi4 DFT -> HOMO/LUMO
└── homo_lumo_results.csv   # appended results (one row per compound/environment run)
bde/                  # O–H bond dissociation energy at every hydroxyl site
figures/              # Summary figures/tables
```

**Note on MD analysis files:** whether the original `.xvg` analysis output
still exists (it may be recoverable from Google Drive) is being checked
separately; `md/` currently holds only the setup files (mdp, topology).
The values in Table 1 below are taken directly from the completed
analysis regardless.

## Method

### 1. Molecular dynamics (non-reactive; proximity/solvation only)

10 ns × 2 independent replicates (independent initial velocities) per
system, for rosavin, rosin, and L-arabinose, each with H₂O₂ present in
solution, using [GROMACS](https://www.gromacs.org/). For each system, the
following were computed:

- Center-of-mass distance to H₂O₂ (nm)
- Number of H-bonds between solute and H₂O₂
- Number of H-bonds between solute and H₂O
- Total solvent-accessible surface area, SASA (nm²)
- Cinnamyl-group SASA (nm²)
- Hydroxyl-group SASA (nm²)

**Note:** MD here is explicitly non-reactive — it captures proximity and
solvation behavior only, not chemical (redox) reactivity.

### 2. Frontier molecular orbitals (DFT)

HOMO/LUMO energies computed with [Psi4](https://psicode.org/) (v1.11) at
the B3LYP/6-31G* level, for rosavin, rosin, and L-arabinose, in both gas
phase and implicit water (PCM solvation). 3D geometries were generated
from SMILES with [RDKit](https://www.rdkit.org/) (MMFF-optimized); this
MMFF geometry was used directly for a **single-point** DFT calculation
(no DFT-level re-optimization).

### 3. O–H bond dissociation energy (BDE)

BDE calculated by DFT (B3LYP/6-31G*, implicit water) at **every** hydroxyl
site — 6 sites in rosavin, 4 in rosin (an exhaustive, unbiased comparison,
not a subset):

```
BDE = E(radical) + E(H atom) − E(parent molecule)
```

where the radical results from homolytic O–H bond cleavage at one
hydroxyl site, with the parent and radical geometries fully optimized.

**Note:** confirmed to use Psi4 (not Gaussian) for the HOMO/LUMO
calculations above; BDE was very likely computed the same way, as part of
the same workflow, but this has not yet been separately confirmed against
the actual notebook code.

## Results

### Table 1: Molecular dynamics (10 ns × 2 replicates per compound)

| Solute | COM distance to H₂O₂ (nm) | H-bonds, solute–H₂O₂ | Total SASA (nm²) | Cinnamyl SASA (nm²) | Hydroxyl SASA (nm²) |
|---|---|---|---|---|---|
| Rosavin | 1.714 ± 0.138 | 0.016 ± 0.008 | 6.622 ± 0.165 | 2.035 ± 0.078 | 2.398 ± 0.104 |
| Rosin | 1.808 ± 0.015 | 0.008 ± 0.005 | 5.539 ± 0.006 | 2.308 ± 0.007 | 1.748 ± 0.006 |
| L-arabinose | 1.499 ± 0.001 | 0.010 ± 0.002 | 3.056 ± 0.001 | — | 1.850 ± 0.000 |

No significant compound-specific difference in direct physical
interaction with H₂O₂ (distance, H-bonding, SASA all track molecular
size rather than reactivity).

### HOMO/LUMO (DFT, B3LYP/6-31G*, Psi4)

| Solute | Environment | HOMO (eV) | LUMO (eV) | Gap (eV) |
|---|---|---|---|---|
| Rosavin | gas | -6.224 | -0.984 | 5.240 |
| Rosin | gas | -6.113 | -0.657 | 5.456 |
| L-arabinose | gas | -6.972 | 1.028 | 8.000 |
| Rosavin | water (PCM) | -6.149 | -0.888 | 5.261 |
| Rosin | water (PCM) | -6.214 | -0.746 | 5.468 |
| L-arabinose | water (PCM) | -7.000 | 1.690 | 8.690 |

The cinnamyl (π-conjugated) group, present in rosavin and rosin but
absent in L-arabinose, confers substantially greater electron-donating
capacity (HOMO closer to zero): L-arabinose's HOMO is 0.75–0.86 eV lower
than rosavin's or rosin's, in both gas phase and water — a difference
exceeding typical DFT uncertainty. Rosavin vs. rosin is a smaller effect
whose direction depends on environment (rosin's HOMO is higher in gas
phase; rosavin's is higher in water).

### O–H bond dissociation energy (BDE)

Within the shared cinnamyl-glycoside scaffold, rosavin has more low-BDE
hydroxyl sites than rosin (4 of 6 vs. 0 of 4 below 110 kcal/mol), and its
most reactive site is ~3.9 kcal/mol lower than rosin's most reactive
site.

**Taken together**, these results point to the cinnamyl moiety as the
dominant electronic driver of oxidizability, with rosavin's greater
number of low-BDE hydroxyl sites as an additional, descriptive
(non-statistically-tested) structural trend distinguishing it from rosin.

## Requirements

- [GROMACS](https://www.gromacs.org/)
- [Psi4](https://psicode.org/)
- [RDKit](https://www.rdkit.org/) (SMILES → 3D geometry, MMFF pre-optimization)

## Citation

Brink DF, Sapp TL, Ghafoor TS, Boyland PA, Tamazawa YC, Kaur G, Shults NV,
Sullivan RD, Suzuki YJ. Antioxidant Properties of a *Rhodiola rosea*
Constituent, Rosavin. *Life* (2026), 16(9):1510.
DOI: [10.3390/life16091510](https://doi.org/10.3390/life16091510)

## License

[MIT License](LICENSE)
