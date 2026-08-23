from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
CSV_PATH = REPO_ROOT / "postProcessing" / "final_kpis.csv"
OUTPUT_DIR = REPO_ROOT / "figures" / "biological_response"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(CSV_PATH)

INLET_DO = 0.50
DO_SAT_REF = 10.8

labels = {
    "A": "Distributed (A)",
    "B": "Central (B)",
    "C": "Lateral (C)",
}

cases = [labels[c] for c in df["Case"]]

x = np.arange(len(cases))
width = 0.36

fig, ax = plt.subplots(figsize=(7, 4.8))

bars_in = ax.bar(
    x - width / 2,
    [INLET_DO] * len(df),
    width,
    label="Inlet",
)

bars_out = ax.bar(
    x + width / 2,
    df["outlet_DO_mgL"],
    width,
    label="Outlet",
)

ax.set_xlabel("Diffuser layout")
ax.set_ylabel("Dissolved oxygen (mg/L)")
ax.set_title("Inlet and Flow-Weighted Outlet Dissolved Oxygen")

ax.set_xticks(x)
ax.set_xticklabels(cases)

ax.bar_label(
    bars_in,
    labels=[f"{INLET_DO:.2f}"] * len(df),
    padding=3,
)

ax.bar_label(
    bars_out,
    labels=[f"{v:.2f}" for v in df["outlet_DO_mgL"]],
    padding=3,
)

ax.axhline(
    DO_SAT_REF,
    linestyle="--",
    linewidth=1.2,
)

ax.set_ylim(
    0,
    max(max(df["outlet_DO_mgL"]), DO_SAT_REF) * 1.12,
)

ax.legend(frameon=False, loc="upper left")

fig.text(
    0.5,
    0.01,
    r"Dashed line: reference surface $DO_{sat}$ at 12 °C ≈ 10.8 mg/L",
    ha="center",
    fontsize=9,
)

fig.tight_layout(rect=[0, 0.06, 1, 1])

fig.savefig(
    OUTPUT_DIR / "05_DO_inlet_outlet_comparison.png",
    dpi=300,
    bbox_inches="tight",
)

plt.show()
