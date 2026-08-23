from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
CSV_PATH = REPO_ROOT / "postProcessing" / "final_kpis.csv"
OUTPUT_DIR = REPO_ROOT / "figures" / "biological_response"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(CSV_PATH)

labels = {
    "A": "Distributed (A)",
    "B": "Central (B)",
    "C": "Lateral (C)",
}

cases = [labels[c] for c in df["Case"]]
outlet_cod = df["outlet_biodegradable_COD_mgL"]
removal = df["biodegradable_COD_removal_pct"]

influent_cod = 271.82  # mg/L, SS + XS

fig, ax = plt.subplots(figsize=(6.5, 4.5))

bars = ax.bar(cases, outlet_cod)

ax.set_xlabel("Diffuser layout")
ax.set_ylabel("Outlet biodegradable COD (mg/L)")
ax.set_title("Flow-Weighted Outlet Biodegradable Substrate COD")

ax.text(
    0.5,
    0.94,
    f"Influent biodegradable substrate COD = {influent_cod:.2f} mg/L",
    transform=ax.transAxes,
    ha="center",
    va="top",
    fontsize=10,
)

for bar, cod, removal_pct in zip(
    bars,
    outlet_cod,
    removal,
):
    ax.text(
        bar.get_x() + bar.get_width() / 2,
        bar.get_height(),
        f"{cod:.2f} mg/L\n{removal_pct:.2f}% reduction",
        ha="center",
        va="bottom",
        fontsize=10,
    )

ax.set_ylim(0, max(outlet_cod) * 1.35)

fig.tight_layout()

fig.savefig(
    OUTPUT_DIR / "03_outlet_biodegradable_COD.png",
    dpi=300,
    bbox_inches="tight",
)

plt.show()
