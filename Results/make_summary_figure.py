import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

compounds = ["Rosavin", "Rosin", "L-arabinose"]
colors = ["#4C72B0", "#DD8452", "#55A868"]

# ---- Table 1: MD (10 ns x 2 replicates) --------------------------------
distance = {"mean": [1.714, 1.808, 1.499], "sd": [0.138, 0.015, 0.001]}
hbond = {"mean": [0.016, 0.008, 0.010], "sd": [0.008, 0.005, 0.002]}
sasa_total = {"mean": [6.622, 5.539, 3.056], "sd": [0.165, 0.006, 0.001]}

# ---- Table 2: HOMO/LUMO (DFT, B3LYP/6-31G*, Psi4) ----------------------
homo_gas = [-6.224, -6.113, -6.972]
homo_water = [-6.149, -6.214, -7.000]
lumo_gas = [-0.984, -0.657, 1.028]
lumo_water = [-0.888, -0.746, 1.690]

fig, axes = plt.subplots(2, 3, figsize=(13, 8))
fig.suptitle("Rosavin vs. Rosin vs. L-arabinose: MD and DFT summary", fontsize=14, fontweight="bold")

# --- Row 1: MD ---
ax = axes[0, 0]
ax.bar(compounds, distance["mean"], yerr=distance["sd"], color=colors, capsize=4)
ax.set_ylabel("nm")
ax.set_title("COM distance to H$_2$O$_2$")

ax = axes[0, 1]
ax.bar(compounds, hbond["mean"], yerr=hbond["sd"], color=colors, capsize=4)
ax.set_title("H-bonds, solute\N{EN DASH}H$_2$O$_2$")

ax = axes[0, 2]
ax.bar(compounds, sasa_total["mean"], yerr=sasa_total["sd"], color=colors, capsize=4)
ax.set_ylabel("nm$^2$")
ax.set_title("Total SASA")

# --- Row 2: DFT ---
x = np.arange(len(compounds))
width = 0.35

ax = axes[1, 0]
ax.bar(x - width/2, homo_gas, width, label="gas", color="#8C8C8C")
ax.bar(x + width/2, homo_water, width, label="water (PCM)", color="#4C72B0")
ax.set_xticks(x); ax.set_xticklabels(compounds)
ax.set_ylabel("eV")
ax.set_title("HOMO")
ax.legend(fontsize=8)
ax.axhline(0, color="black", linewidth=0.6)

ax = axes[1, 1]
ax.bar(x - width/2, lumo_gas, width, label="gas", color="#8C8C8C")
ax.bar(x + width/2, lumo_water, width, label="water (PCM)", color="#4C72B0")
ax.set_xticks(x); ax.set_xticklabels(compounds)
ax.set_title("LUMO")
ax.legend(fontsize=8)
ax.axhline(0, color="black", linewidth=0.6)

ax = axes[1, 2]
gap_gas = [l - h for l, h in zip(lumo_gas, homo_gas)]
gap_water = [l - h for l, h in zip(lumo_water, homo_water)]
ax.bar(x - width/2, gap_gas, width, label="gas", color="#8C8C8C")
ax.bar(x + width/2, gap_water, width, label="water (PCM)", color="#4C72B0")
ax.set_xticks(x); ax.set_xticklabels(compounds)
ax.set_title("HOMO\N{EN DASH}LUMO gap")
ax.legend(fontsize=8)

for row in axes:
    for ax in row:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(axis="x", labelrotation=15)

fig.text(0.5, 0.01,
         "BDE (O\u2013H bond dissociation energy) not shown \u2014 only 2 of 20 site/environment\n"
         "values are confirmed against source notebooks so far.",
         ha="center", fontsize=9, style="italic", color="#555555")

plt.tight_layout(rect=[0, 0.04, 1, 0.96])
plt.savefig("results_summary.png", dpi=200, bbox_inches="tight")
print("saved results_summary.png")
