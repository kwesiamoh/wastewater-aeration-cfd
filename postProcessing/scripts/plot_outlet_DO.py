from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
CSV_PATH = REPO_ROOT / "postProcessing" / "final_kpis.csv"
OUTPUT_DIR = REPO_ROOT / "figures" / "biological_response"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(CSV_PATH)

DO_SAT_REF = 10.8  # mg/L, reference surface saturation at 12 °C

labels = {
    "A": "Distributed (A)",
    "B": "Central (B)",
    "C": "Lateral (C)",
}

cases = [labels[c] for c in df["Case"]]
outlet_do = df["outlet_DO_mgL"]

fig, ax = plt.subplots(figsize=(6.5, 4.5))

bars = ax.bar(cases, outlet_do)

ax.axhline(
    DO_SAT_REF,
    linestyle="--",
    linewidth=1.2,
)

ax.set_xlabel("Diffuser layout")
ax.set_ylabel("Outlet dissolved oxygen (mg/L)")
ax.set_title("Flow-Weighted Outlet Dissolved Oxygen")

ax.bar_label(
    bars,
    labels=[f"{v:.2f}" for v in outlet_do],
    padding=3,
)

ax.set_ylim(0, max(max(outlet_do), DO_SAT_REF) * 1.12)

fig.text(
    0.5,
    0.01,
    r"Dashed line: reference surface $DO_{sat}$ at 12 °C ≈ 10.8 mg/L",
    ha="center",
    fontsize=9,
)

fig.tight_layout(rect=[0, 0.06, 1, 1])

fig.savefig(
    OUTPUT_DIR / "04_outlet_DO_comparison.png",
    dpi=300,
    bbox_inches="tight",
)

plt.show()
