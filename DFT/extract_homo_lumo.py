"""
Extract HOMO/LUMO orbital energies from Gaussian .log output.

Requires the route line to include `pop=full` (or `pop=reg`), which prints
"Alpha occ. eigenvalues" / "Alpha virt. eigenvalues" blocks. The HOMO is
the last occupied eigenvalue printed; the LUMO is the first virtual
eigenvalue printed, from the FINAL such block in the file (i.e. from the
optimized geometry, not an intermediate SCF step).

Usage:
    python extract_homo_lumo.py gas_phase/rosavin.log gas_phase/rosin.log ...
"""
import sys
import re
from pathlib import Path

HARTREE_TO_EV = 27.2114

OCC_RE = re.compile(r"Alpha\s+occ\.\s+eigenvalues\s+--\s+(.*)")
VIRT_RE = re.compile(r"Alpha\s+virt\.\s+eigenvalues\s+--\s+(.*)")


def parse_homo_lumo(log_path):
    """Return (HOMO_eV, LUMO_eV, gap_eV) from the LAST orbital block in the file."""
    occ_energies = []
    virt_energies = []
    last_occ, last_virt = None, None

    with open(log_path) as fh:
        for line in fh:
            m_occ = OCC_RE.search(line)
            m_virt = VIRT_RE.search(line)
            if m_occ:
                occ_energies.extend(float(x) for x in m_occ.group(1).split())
            elif m_virt:
                virt_energies.extend(float(x) for x in m_virt.group(1).split())
            elif occ_energies or virt_energies:
                # a non-eigenvalue line means this orbital block just ended
                if occ_energies and virt_energies:
                    last_occ, last_virt = occ_energies[-1], virt_energies[0]
                occ_energies, virt_energies = [], []

    if occ_energies and virt_energies:  # in case the file ends mid-block
        last_occ, last_virt = occ_energies[-1], virt_energies[0]

    if last_occ is None or last_virt is None:
        raise ValueError(
            f"No complete HOMO/LUMO block found in {log_path}. "
            "Check that the route line includes 'pop=full' and the job finished."
        )

    homo_ev = last_occ * HARTREE_TO_EV
    lumo_ev = last_virt * HARTREE_TO_EV
    return homo_ev, lumo_ev, lumo_ev - homo_ev


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    print(f"{'File':<40}{'HOMO (eV)':>12}{'LUMO (eV)':>12}{'Gap (eV)':>12}")
    print("-" * 76)
    for path in sys.argv[1:]:
        try:
            homo, lumo, gap = parse_homo_lumo(path)
            print(f"{Path(path).name:<40}{homo:>12.3f}{lumo:>12.3f}{gap:>12.3f}")
        except ValueError as e:
            print(f"{Path(path).name:<40}  ERROR: {e}")
