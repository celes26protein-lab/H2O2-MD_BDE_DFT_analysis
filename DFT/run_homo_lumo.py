"""
Rosavin / rosin / L-arabinose -- HOMO/LUMO via Psi4 (B3LYP/6-31G*)

Builds a 3D geometry from a SMILES string with RDKit (MMFF-optimized as a
starting guess), then runs a Psi4 DFT geometry optimization + single-point
to get HOMO/LUMO orbital energies, in either gas phase or implicit water
(PCM). Appends one row per run to homo_lumo_results.csv.

This version's Psi4 logic (single-point, no re-optimization; orbital
extraction via epsilon_a_subset("AO","ALL") indexed by nalpha()) mirrors
a confirmed, already-run version of this calculation for rosavin (gas
phase), which reproduced HOMO = -6.224 eV, LUMO = -0.984 eV, gap = 5.240
eV. It has not been re-run in this environment (psi4 is conda-only and
was not available here to test).

Usage:
    python run_homo_lumo.py <name> "<SMILES>" <gas|water>

Example:
    python run_homo_lumo.py rosavin \
        'O[C@H]1CO[C@@H](OC[C@H]2O[C@@H](OC\\C=C\\C3=CC=CC=C3)[C@H](O)[C@@H](O)[C@@H]2O)[C@H](O)[C@H]1O' \
        gas
"""
import sys
import csv
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

import psi4

HARTREE_TO_EV = 27.2114
RESULTS_CSV = Path(__file__).parent / "homo_lumo_results.csv"


def smiles_to_xyz_block(smiles, seed=42):
    """SMILES -> 3D geometry (RDKit/MMFF), atom lines only (no charge/mult line)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles}")
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, randomSeed=seed) != 0:
        raise RuntimeError("RDKit 3D embedding failed")
    AllChem.MMFFOptimizeMolecule(mol)

    conf = mol.GetConformer()
    lines = []
    for atom in mol.GetAtoms():
        p = conf.GetAtomPosition(atom.GetIdx())
        lines.append(f"{atom.GetSymbol()} {p.x:.6f} {p.y:.6f} {p.z:.6f}")
    return "\n".join(lines), mol.GetNumAtoms(), Chem.rdMolDescriptors.CalcMolFormula(mol)


def run_homo_lumo(name, smiles, environment, memory="4 GB", nthreads=2):
    assert environment in ("gas", "water")

    xyz_block, n_atoms, formula = smiles_to_xyz_block(smiles)
    print(f"{name}: {n_atoms} atoms, formula {formula}")

    psi4.core.set_output_file(f"{name}_{environment}.log", False)
    psi4.set_memory(memory)
    psi4.set_num_threads(nthreads)

    psi4_geom = f"""
0 1
{xyz_block}
units angstrom
"""
    molecule = psi4.geometry(psi4_geom)

    options = {"basis": "6-31G*", "reference": "rhf"}
    if environment == "water":
        # Implicit solvent (PCM). Requires Psi4 built with PCMSolver.
        options.update({
            "pcm": True,
            "pcm_scf_type": "total",
        })
        psi4.pcm_helper("""
            Units = Angstrom
            Medium {
                SolverType = IEFPCM
                Solvent = Water
            }
            Cavity {
                RadiiSet = UFF
                Type = GePol
                Scaling = False
                Area = 0.3
                Mode = Implicit
            }
        """)
    psi4.set_options(options)

    # Single-point only -- the RDKit/MMFF geometry is used as-is, with no
    # DFT-level re-optimization (matches the confirmed working method).
    energy, wfn = psi4.energy("b3lyp", return_wfn=True, molecule=molecule)

    homo_idx = wfn.nalpha() - 1   # 0-indexed
    lumo_idx = wfn.nalpha()
    eps = wfn.epsilon_a_subset("AO", "ALL").to_array()
    homo_ev = eps[homo_idx] * HARTREE_TO_EV
    lumo_ev = eps[lumo_idx] * HARTREE_TO_EV
    gap_ev = lumo_ev - homo_ev

    print(f"  Total Energy: {energy:.6f} Hartree")
    print(f"  HOMO = {homo_ev:.3f} eV   LUMO = {lumo_ev:.3f} eV   gap = {gap_ev:.3f} eV")

    write_header = not RESULTS_CSV.exists()
    with open(RESULTS_CSV, "a", newline="") as fh:
        writer = csv.writer(fh)
        if write_header:
            writer.writerow(["compound", "environment", "n_atoms", "formula",
                              "HOMO_eV", "LUMO_eV", "gap_eV"])
        writer.writerow([name, environment, n_atoms, formula,
                          f"{homo_ev:.4f}", f"{lumo_ev:.4f}", f"{gap_ev:.4f}"])

    return homo_ev, lumo_ev, gap_ev


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(__doc__)
        sys.exit(1)
    name, smiles, environment = sys.argv[1], sys.argv[2], sys.argv[3]
    run_homo_lumo(name, smiles, environment)
