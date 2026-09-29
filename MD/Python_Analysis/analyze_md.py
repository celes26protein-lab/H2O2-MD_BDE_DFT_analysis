"""
Rosavin / rosin / L-arabinose vs H2O2 -- MD analysis (Table 1)

Reads GROMACS .xvg output (gmx distance, gmx hbond, gmx sasa) for each
compound and each of the two independent 10-ns replicates, and reports
mean +/- SD for:

    - Center-of-mass distance to H2O2 (nm)
    - Number of H-bonds between solute and H2O2
    - Total SASA (nm^2)
    - Cinnamyl-group SASA (nm^2)      [not applicable to L-arabinose]
    - Hydroxyl-group SASA (nm^2)

Expected input files, in the same directory as this script (xvg/):

    distance_h2o2_<compound>_rep<1|2>.xvg
    hbond_h2o2_<compound>_rep<1|2>.xvg
    sasa_total_<compound>_rep<1|2>.xvg
    sasa_cinnamyl_<compound>_rep<1|2>.xvg     (rosavin, rosin only)
    sasa_hydroxyl_<compound>_rep<1|2>.xvg

where <compound> is one of: rosavin, rosin, l_arabinose

NOTE on averaging convention: the two replicates are pooled (all frames
from rep1 and rep2 combined) before computing mean +/- SD, matching
"calculated across two independent 10 ns production simulations" in the
paper's table footnote. If your original numbers were instead computed as
the mean/SD *across the two replicate means* (n=2), change
POOL_FRAMES to False below.
"""
import numpy as np
from pathlib import Path

XVG_DIR = Path(__file__).parent / "xvg"
COMPOUNDS = ["rosavin", "rosin", "l_arabinose"]
REPLICATES = [1, 2]
POOL_FRAMES = True   # see note above


def read_xvg(path):
    """Read the last numeric column of a GROMACS .xvg file (skips # and @ lines)."""
    values = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(("#", "@")):
                continue
            cols = line.split()
            if cols:
                values.append(float(cols[-1]))
    return np.array(values)


def load_metric(prefix, compound):
    """Load and combine both replicates for one metric/compound. Returns None if missing."""
    arrays = []
    for rep in REPLICATES:
        f = XVG_DIR / f"{prefix}_{compound}_rep{rep}.xvg"
        if not f.exists():
            return None
        arrays.append(read_xvg(f))

    if POOL_FRAMES:
        return np.concatenate(arrays)
    else:
        means = np.array([a.mean() for a in arrays])
        return means  # mean/std of this will be computed from n=2 replicate means


def summarize(prefix, compound, label, unit=""):
    data = load_metric(prefix, compound)
    if data is None:
        return None
    return data.mean(), data.std()


METRICS = [
    ("distance_h2o2",  "COM distance to H2O2",        "nm"),
    ("hbond_h2o2",     "H-bonds, solute-H2O2",         ""),
    ("sasa_total",     "Total SASA",                   "nm^2"),
    ("sasa_cinnamyl",  "Cinnamyl-group SASA",          "nm^2"),
    ("sasa_hydroxyl",  "Hydroxyl-group SASA",          "nm^2"),
]

print(f"{'Compound':<14}{'Metric':<26}{'Mean':>10}{'SD':>10}   Unit")
print("-" * 68)

results = {}
for compound in COMPOUNDS:
    results[compound] = {}
    for prefix, label, unit in METRICS:
        r = summarize(prefix, compound, label, unit)
        results[compound][prefix] = r
        if r is None:
            print(f"{compound:<14}{label:<26}{'(not found)':>10}")
        else:
            mean, sd = r
            print(f"{compound:<14}{label:<26}{mean:>10.3f}{sd:>10.3f}   {unit}")
    print()

# ---- Table 1-style summary --------------------------------------------
print("\n=== Table 1 (mean +/- SD) ===\n")
header = f"{'Solute':<14}" + "".join(f"{label:>22}" for _, label, _ in METRICS)
print(header)
for compound in COMPOUNDS:
    row = f"{compound:<14}"
    for prefix, _, _ in METRICS:
        r = results[compound][prefix]
        cell = f"{r[0]:.3f} +/- {r[1]:.3f}" if r is not None else "--"
        row += f"{cell:>22}"
    print(row)
