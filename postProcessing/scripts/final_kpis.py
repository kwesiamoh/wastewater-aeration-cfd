from pathlib import Path
import pandas as pd


# ---------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]

OUTPUT_CSV = REPO_ROOT / "postProcessing" / "final_kpis.csv"

AERATION_CASES = {
    "A": REPO_ROOT / "cases" / "aeration" / "CaseA_Distributed",
    "B": REPO_ROOT / "cases" / "aeration" / "CaseB_Central",
    "C": REPO_ROOT / "cases" / "aeration" / "CaseC_Lateral",
}

BIOLOGICAL_CASES = {
    "A": REPO_ROOT / "cases" / "biological" / "CaseA",
    "B": REPO_ROOT / "cases" / "biological" / "CaseB",
    "C": REPO_ROOT / "cases" / "biological" / "CaseC",
}

# Use one consistent developed-flow averaging window
AERATION_START = 180.0
AERATION_END = 300.0


# ---------------------------------------------------------------------
# OpenFOAM postProcessing readers
# ---------------------------------------------------------------------

def read_openfoam_table(path):
    """
    Read a simple OpenFOAM postProcessing .dat file.

    Returns
    -------
    columns : list[str]
    rows    : list[list[float]]
    """
    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(f"Missing postProcessing file:\n{path}")

    columns = None
    rows = []

    with path.open() as f:
        for line in f:
            stripped = line.strip()

            if not stripped:
                continue

            if stripped.startswith("# Time"):
                columns = stripped.lstrip("#").split()
                continue

            if stripped.startswith("#"):
                continue

            parts = stripped.split()

            try:
                rows.append([float(x) for x in parts])
            except ValueError:
                continue

    if columns is None:
        raise ValueError(f"Could not find '# Time' header in:\n{path}")

    if not rows:
        raise ValueError(f"No numeric data found in:\n{path}")

    return columns, rows


def get_column(path, column_name):
    columns, rows = read_openfoam_table(path)

    try:
        idx = columns.index(column_name)
    except ValueError as exc:
        raise ValueError(
            f"Column '{column_name}' not found in {path}\n"
            f"Available columns: {columns}"
        ) from exc

    return [(row[0], row[idx]) for row in rows]


def mean_over_window(path, column_name, start, end):
    values = [
        value
        for time, value in get_column(path, column_name)
        if start <= time <= end
    ]

    if not values:
        raise ValueError(
            f"No values for '{column_name}' between "
            f"{start} and {end} in {path}"
        )

    return sum(values) / len(values)


def latest_value(path, column_name):
    values = get_column(path, column_name)
    return values[-1][1]


# ---------------------------------------------------------------------
# Aeration KPIs
# ---------------------------------------------------------------------

def extract_aeration(case_path):
    pp = case_path / "postProcessing"

    kla_s_inv = mean_over_window(
        pp / "tankKLa" / "0" / "volFieldValue.dat",
        "weightedVolAverage(kLa)",
        AERATION_START,
        AERATION_END,
    )

    gas_fraction = mean_over_window(
        pp / "tankAir" / "0" / "volFieldValue.dat",
        "weightedVolAverage(alpha.air)",
        AERATION_START,
        AERATION_END,
    )

    low_velocity_fraction = mean_over_window(
        pp / "lowVelocityFraction" / "0" / "volFieldValue.dat",
        "weightedVolAverage(lowVelocityMask)",
        AERATION_START,
        AERATION_END,
    )

    return {
        # kLa output is s^-1; convert to h^-1
        "kLa_h-1": kla_s_inv * 3600.0,

        # Fractions -> percent
        "gas_holdup_pct": gas_fraction * 100.0,
        "low_velocity_pct": low_velocity_fraction * 100.0,
    }


# ---------------------------------------------------------------------
# Reduced ASM1 KPIs
# ---------------------------------------------------------------------

def extract_biological(case_path):
    outlet_file = (
        case_path
        / "postProcessing"
        / "outletFlowWeighted"
        / "300"
        / "surfaceFieldValue.dat"
    )

    # ASM1 concentration fields are kg/m3.
    # 1 kg/m3 = 1000 mg/L.
    to_mgL = 1000.0

    ss = latest_value(
        outlet_file,
        "weightedAverage(SS)"
    ) * to_mgL

    xs = latest_value(
        outlet_file,
        "weightedAverage(XS)"
    ) * to_mgL

    do = latest_value(
        outlet_file,
        "weightedAverage(SO)"
    ) * to_mgL

    substrate_cod = latest_value(
        outlet_file,
        "weightedAverage(substrateCOD)"
    ) * to_mgL

    return {
        "outlet_DO_mgL": do,
        "outlet_Ss_mgL": ss,
        "outlet_Xs_mgL": xs,
        "outlet_biodegradable_COD_mgL": substrate_cod,
    }


def get_influent_substrate_cod():
    """
    The three biological cases use the same influent composition.
    Read it directly from the final Case A reduced-ASM1 inlet result.
    """
    inlet_file = (
        BIOLOGICAL_CASES["A"]
        / "postProcessing"
        / "inletASM"
        / "300"
        / "surfaceFieldValue.dat"
    )

    return (
        latest_value(
            inlet_file,
            "areaAverage(substrateCOD)"
        )
        * 1000.0
    )


# ---------------------------------------------------------------------
# Build final KPI table
# ---------------------------------------------------------------------

influent_cod = get_influent_substrate_cod()

records = []

for case in ["A", "B", "C"]:
    aeration = extract_aeration(AERATION_CASES[case])
    biological = extract_biological(BIOLOGICAL_CASES[case])

    outlet_cod = biological["outlet_biodegradable_COD_mgL"]

    removal = (
        (influent_cod - outlet_cod)
        / influent_cod
        * 100.0
    )

    records.append(
        {
            "Case": case,
            **aeration,
            **biological,
            "influent_biodegradable_COD_mgL": influent_cod,
            "biodegradable_COD_removal_pct": removal,
        }
    )


df = pd.DataFrame(records)

df.to_csv(OUTPUT_CSV, index=False)

print()
print(
    f"Aeration averaging window: "
    f"{AERATION_START:.0f}-{AERATION_END:.0f} s"
)

print(
    f"Influent biodegradable COD: "
    f"{influent_cod:.5f} mg/L"
)

print()
print(df.to_string(index=False))

print()
print(f"Saved: {OUTPUT_CSV}")
