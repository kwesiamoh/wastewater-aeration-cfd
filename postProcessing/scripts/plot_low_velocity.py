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
low_velocity = df["low_velocity_pct"]

fig, ax = plt.subplots(figsize=(10, 6))

bars = ax.bar(cases, low_velocity, width=0.58, zorder=3)

ax.set_title("Low-velocity volume fraction", fontsize=16, pad=12)
ax.set_xlabel("Diffuser layout", fontsize=13, labelpad=8)
ax.set_ylabel("Low-velocity volume fraction (%)", fontsize=13, labelpad=8)

ax.tick_params(axis="both", labelsize=11)
ax.set_ylim(0, max(low_velocity) * 1.15)
ax.yaxis.grid(True, linewidth=0.8, zorder=0)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

for bar, value in zip(bars, low_velocity):
    ax.text(
        bar.get_x() + bar.get_width() / 2,
        value + max(low_velocity) * 0.02,
        f"{value:.1f}%",
        ha="center",
        va="bottom",
        fontsize=11,
    )

fig.text(
    0.5,
    0.01,
    r"Low velocity: $|U_{water}| < 0.2$ m s$^{-1}$; developed period: 180–300 s",
    ha="center",
    fontsize=9,
)

plt.tight_layout(rect=[0, 0.05, 1, 1])

output = SCRIPT_DIR.parent.parent / "figures" / "hydrodynamics"
output.mkdir(parents=True, exist_ok=True)

plt.savefig(
    output / "02_low_velocity_comparison.png",
    dpi=300,
    bbox_inches="tight",
)

plt.show()
