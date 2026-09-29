# DFT: HOMO/LUMO (Frontier Molecular Orbitals)

Computes HOMO/LUMO orbital energies for rosavin, rosin, and L-arabinose,
in gas phase and implicit water, to compare their intrinsic
electron-donating capacity (oxidizability) as part of the
[Rosavin–H₂O₂ mechanism study](../README.md).

## Method

1. **3D geometry:** built from a SMILES string with [RDKit](https://www.rdkit.org/)
   (`EmbedMolecule` + `MMFFOptimizeMolecule`). This MMFF-optimized geometry
   is used directly — there is **no DFT-level re-optimization**.
2. **Single-point DFT:** [Psi4](https://psicode.org/) (v1.11), B3LYP/6-31G*,
   RHF reference.
3. **Implicit water (PCM):** IEFPCM solvent model via `psi4.pcm_helper`
   (UFF radii, GePol cavity, `Scaling = False`, `Area = 0.3`).
4. **HOMO/LUMO extraction:** `wfn.epsilon_a_subset("AO", "ALL")`, indexed
   at `nalpha() - 1` (HOMO) and `nalpha()` (LUMO). Hartree → eV via
   ×27.2114.

## Contents

```
dft_homo_lumo/
├── DFT_README.md            # this file
├── run_homo_lumo.py         # SMILES -> RDKit 3D geometry -> Psi4 -> HOMO/LUMO
└── *.log                    # Psi4 output logs (one per compound/environment run)
```

**Not included:** `homo_lumo_results.csv` — omitted due to file size,
along with the equivalent BDE results CSV (see
[../bde/BDE_README.md](../bde/BDE_README.md)). The completed HOMO/LUMO
values are reported below and in the main README instead.

## Usage

```bash
python3 run_homo_lumo.py <name> "<SMILES>" <gas|water>
```

Run once per compound, per environment (6 runs total):

```bash
python3 run_homo_lumo.py rosavin      "O[C@H]1CO[C@@H](OC[C@H]2O[C@@H](OC\C=C\C3=CC=CC=C3)[C@H](O)[C@@H](O)[C@@H]2O)[C@H](O)[C@H]1O" gas
python3 run_homo_lumo.py rosavin      "O[C@H]1CO[C@@H](OC[C@H]2O[C@@H](OC\C=C\C3=CC=CC=C3)[C@H](O)[C@@H](O)[C@@H]2O)[C@H](O)[C@H]1O" water

python3 run_homo_lumo.py rosin        "C1=CC=C(C=C1)/C=C/CO[C@H]2[C@@H]([C@H]([C@@H]([C@H](O2)CO)O)O)O" gas
python3 run_homo_lumo.py rosin        "C1=CC=C(C=C1)/C=C/CO[C@H]2[C@@H]([C@H]([C@@H]([C@H](O2)CO)O)O)O" water

python3 run_homo_lumo.py l_arabinose  "C1[C@@H]([C@@H]([C@H](C(O1)O)O)O)O" gas
python3 run_homo_lumo.py l_arabinose  "C1[C@@H]([C@@H]([C@H](C(O1)O)O)O)O" water
```

Each run appends one row to `homo_lumo_results.csv` (compound,
environment, atom count, formula, HOMO, LUMO, gap).

## Requirements

- [Psi4](https://psicode.org/) (conda-only: `conda install -c conda-forge -c psi4 psi4`)
- [RDKit](https://www.rdkit.org/)

## Results

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

## Caveats

- **No DFT-level geometry optimization.** HOMO/LUMO can be sensitive to
  geometry; using the MMFF structure directly is a practical
  simplification, not a validated equivalence. If revisiting, re-running
  at least rosavin with a B3LYP-optimized geometry and confirming the
  cinnamyl-vs-arabinose gap holds would strengthen this result.
- Single conformer per compound (no conformational sampling).
