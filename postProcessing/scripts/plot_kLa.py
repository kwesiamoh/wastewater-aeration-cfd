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
kla = df["kLa_h-1"]

fig, ax = plt.subplots(figsize=(10, 6))

bars = ax.bar(cases, kla, width=0.58, zorder=3)

ax.set_title(r"Mean oxygen-transfer coefficient, $k_La$", fontsize=16, pad=12)
ax.set_xlabel("Diffuser layout", fontsize=13, labelpad=8)
ax.set_ylabel(r"Mean $k_La$ (h$^{-1}$)", fontsize=13, labelpad=8)

ax.tick_params(axis="both", labelsize=11)
ax.set_ylim(0, max(kla) * 1.18)
ax.yaxis.grid(True, linewidth=0.8, zorder=0)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

for bar, value in zip(bars, kla):
    ax.text(
        bar.get_x() + bar.get_width() / 2,
        value + max(kla) * 0.02,
        f"{value:.2f}",
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

output = SCRIPT_DIR.parent.parent / "figures" / "oxygen_transfer"
output.mkdir(parents=True, exist_ok=True)

plt.savefig(
    output / "01_kLa_comparison.png",
    dpi=300,
    bbox_inches="tight",
)

plt.show()
