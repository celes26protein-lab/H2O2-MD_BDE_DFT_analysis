"""
Rosavin / rosin -- O-H bond dissociation energy (BDE) via Psi4

For every hydroxyl (O-H) site in a compound:
    1. Build the parent molecule from SMILES (RDKit, MMFF-optimized),
       single-point DFT energy (RHF reference, closed shell).
    2. Remove that one H atom, set a radical electron on the O, re-embed
       and MMFF-optimize the resulting doublet radical, then single-point
       DFT energy (UHF reference).
    3. BDE = E(radical) + E(H atom, doublet, UHF) - E(parent)

Matches a confirmed, already-run version of this calculation:
    Rosavin O1-H  (gas):   BDE = 0.172221 Hartree = 108.07 kcal/mol
    Rosavin O23-H (water): BDE = 0.185961 Hartree = 116.69 kcal/mol
This script has NOT been re-run end-to-end in this environment (psi4 is
conda-only and unavailable here); only the RDKit radical-generation logic
has been tested standalone. Please verify against the values above
before trusting results for other sites/compounds.

L-arabinose is not included here -- BDE was reported for rosavin and
rosin only (see main README).

Usage:
    python run_bde.py <name> "<SMILES>" <gas|water>

Example:
    python run_bde.py rosavin \
        'O[C@H]1CO[C@@H](OC[C@H]2O[C@@H](OC\\C=C\\C3=CC=CC=C3)[C@H](O)[C@@H](O)[C@@H]2O)[C@H](O)[C@H]1O' \
        gas
"""
import sys
import csv
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, RWMol

import psi4

HARTREE_TO_KCAL = 627.5095
RESULTS_CSV = Path(__file__).parent / "bde_results.csv"


def build_parent(smiles, seed=42):
    """SMILES -> RDKit mol with a 3D, MMFF-optimized conformer."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles}")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=seed)
    AllChem.MMFFOptimizeMolecule(mol)
    return mol


def find_hydroxyl_sites(mol):
    """Return [(O_idx, H_idx), ...] for every O atom bonded to exactly one H
    (i.e. every -OH oxygen, whether alcohol, phenol, or anomeric)."""
    sites = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "O":
            continue
        h_neighbors = [n.GetIdx() for n in atom.GetNeighbors() if n.GetSymbol() == "H"]
        if len(h_neighbors) == 1:
            sites.append((atom.GetIdx(), h_neighbors[0]))
    return sites


def make_radical_mol(parent_mol, o_idx, h_idx, seed=42):
    """Remove one hydroxyl H, set a radical electron on its O, and
    re-embed + MMFF-optimize a fresh 3D conformer for the resulting
    doublet radical."""
    rw = RWMol(parent_mol)
    rw.RemoveAtom(h_idx)

    # after removing an atom, indices above h_idx shift down by 1
    o_idx_shifted = o_idx - 1 if o_idx > h_idx else o_idx

    rw.GetAtomWithIdx(o_idx_shifted).SetNumRadicalElectrons(1)
    rw.GetAtomWithIdx(o_idx_shifted).SetNoImplicit(True)

    radical = rw.GetMol()
    Chem.SanitizeMol(radical, sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_KEKULIZE)

    AllChem.EmbedMolecule(radical, randomSeed=seed)
    AllChem.MMFFOptimizeMolecule(radical)
    return radical


def mol_to_xyz_block(mol):
    conf = mol.GetConformer()
    lines = []
    for atom in mol.GetAtoms():
        p = conf.GetAtomPosition(atom.GetIdx())
        lines.append(f"{atom.GetSymbol()} {p.x:.6f} {p.y:.6f} {p.z:.6f}")
    return "\n".join(lines)


def set_pcm(environment):
    if environment == "water":
        psi4.set_options({"pcm": True, "pcm_scf_type": "total"})
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


def psi4_single_point(xyz_block, charge, multiplicity, reference, environment, logfile):
    psi4.core.set_output_file(logfile, False)
    reference_upper = reference.upper()
    psi4.set_options({"basis": "6-31G*", "reference": reference_upper})
    set_pcm(environment)
    geom = f"{charge} {multiplicity}\n{xyz_block}\nunits angstrom"
    psi4.geometry(geom)
    return psi4.energy("b3lyp")


def run_bde(name, smiles, environment, memory="4 GB", nthreads=2):
    assert environment in ("gas", "water")
    psi4.set_memory(memory)
    psi4.set_num_threads(nthreads)

    # ---- H atom reference (doublet, UHF) -------------------------------
    h_energy = psi4_single_point(
        "H 0.0 0.0 0.0", charge=0, multiplicity=2, reference="uhf",
        environment=environment, logfile=f"h_atom_{environment}.log",
    )
    print(f"H atom ({environment}): {h_energy:.6f} Hartree")

    # ---- parent molecule (closed shell, RHF) ---------------------------
    parent_mol = build_parent(smiles)
    parent_energy = psi4_single_point(
        mol_to_xyz_block(parent_mol), charge=0, multiplicity=1, reference="rhf",
        environment=environment, logfile=f"{name}_parent_{environment}.log",
    )
    print(f"{name} parent ({environment}): {parent_energy:.6f} Hartree")

    # ---- one radical per hydroxyl site (doublet, UHF) -------------------
    sites = find_hydroxyl_sites(parent_mol)
    print(f"{name}: {len(sites)} hydroxyl site(s) found: {sites}")

    write_header = not RESULTS_CSV.exists()
    with open(RESULTS_CSV, "a", newline="") as fh:
        writer = csv.writer(fh)
        if write_header:
            writer.writerow(["compound", "environment", "O_idx", "H_idx",
                              "parent_hartree", "radical_hartree", "H_atom_hartree",
                              "BDE_hartree", "BDE_kcal_per_mol"])

        for o_idx, h_idx in sites:
            radical_mol = make_radical_mol(parent_mol, o_idx, h_idx)
            radical_energy = psi4_single_point(
                mol_to_xyz_block(radical_mol), charge=0, multiplicity=2, reference="uhf",
                environment=environment, logfile=f"{name}_radical_O{o_idx}_{environment}.log",
            )
            bde_hartree = radical_energy + h_energy - parent_energy
            bde_kcal = bde_hartree * HARTREE_TO_KCAL
            print(f"  O{o_idx}-H{h_idx}: BDE = {bde_hartree:.6f} Hartree = {bde_kcal:.2f} kcal/mol")

            writer.writerow([name, environment, o_idx, h_idx,
                              f"{parent_energy:.6f}", f"{radical_energy:.6f}", f"{h_energy:.6f}",
                              f"{bde_hartree:.6f}", f"{bde_kcal:.2f}"])


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(__doc__)
        sys.exit(1)
    name, smiles, environment = sys.argv[1], sys.argv[2], sys.argv[3]
    run_bde(name, smiles, environment)
