from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

SCRIPT_DIR = Path(__file__).resolve().parent
CSV_PATH = SCRIPT_DIR.parent / "final_kpis.csv"

df = pd.read_csv(CSV_PATH)

labels = {
    "A": "Distributed (A)",
    "B": "Central (B)",
    "C": "Lateral (C)",
}

cases = [labels[c] for c in df["Case"]]
gas_holdup = df["gas_holdup_pct"]

fig, ax = plt.subplots(figsize=(10, 6))

bars = ax.bar(cases, gas_holdup, width=0.58, zorder=3)

ax.set_title("Mean gas holdup", fontsize=16, pad=12)
ax.set_xlabel("Diffuser layout", fontsize=13, labelpad=8)
ax.set_ylabel("Mean gas holdup (%)", fontsize=13, labelpad=8)

ax.tick_params(axis="both", labelsize=11)
ax.set_ylim(0, max(gas_holdup) * 1.18)
ax.yaxis.grid(True, linewidth=0.8, zorder=0)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

for bar, value in zip(bars, gas_holdup):
    ax.text(
        bar.get_x() + bar.get_width() / 2,
        value + max(gas_holdup) * 0.02,
        f"{value:.3f}%",
        ha="center",
        va="bottom",
        fontsize=11,
    )

fig.text(
    0.5,
    0.01,
    "Time-averaged over the developed aeration period (180–300 s)",
    ha="center",
    fontsize=9,
)

plt.tight_layout(rect=[0, 0.05, 1, 1])

output = SCRIPT_DIR.parent.parent / "figures" / "hydrodynamics"
output.mkdir(parents=True, exist_ok=True)

plt.savefig(
    output / "03_gas_holdup_comparison.png",
    dpi=300,
    bbox_inches="tight",
)

plt.show()
