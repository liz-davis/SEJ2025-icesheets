#!/usr/bin/env python3
"""Offline preparation of IMBIE v3 constraints for the SEJ2025 FACTS2 module.

Reproduces original core_imbie_helpers.py conventions:
- common monthly record Jan 1979–Dec 2023
- sign reversal of IMBIE cumulative mass balance anomaly (mm -> SLR cm)
- 1995–2014 monthly-reference mean
- four-year-mean quadratic fit (1979–2022) evaluated at Dec 2023
- observed December endpoint rate via np.gradient(..., edge_order=2)
- SEJ 2005-to-1995–2014 additive reference shifts
- December observations for 2000–2023 (NOT annual means)

No uncertainty sampling is done; this reproduces the deterministic historical
boundary conditions used in the archived SEJ annual interpolation.
"""

from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

COMPONENTS = {
    "E": "imbie3_east_antarctica_mm_partitioned.csv",
    "W": "imbie3_west_antarctica_mm_partitioned.csv",
    "G": "imbie3_greenland_mm_partitioned.csv",
}
CUMULATIVE_COLUMN = "Cumulative mass balance anomaly (mm)"
UNCERTAINTY_COLUMN = "Cumulative mass balance anomaly uncertainty (mm)"
FIT_START_YEAR, FIT_END_YEAR, LAST_OBS_YEAR = 1979, 2022, 2023
REFERENCE_START_YEAR, REFERENCE_END_YEAR = 1995, 2014
MM_TO_CM = 0.1


def decimal_year(dates):
    dates = pd.DatetimeIndex(dates)
    starts = pd.to_datetime(dates.year.astype(str) + "-01-01")
    ends = pd.to_datetime((dates.year + 1).astype(str) + "-01-01")
    return np.asarray(
        dates.year.astype(float)
        + (dates - starts).days.astype(float) / (ends - starts).days.astype(float),
        dtype=float,
    )


def find_header(path):
    with path.open("r", encoding="utf-8-sig", errors="replace") as stream:
        for i, line in enumerate(stream):
            if line.lstrip().startswith("Date,"):
                return i
    raise ValueError(f"Could not locate header starting Date, in {path}")


def read_component(root, component):
    path = root / COMPONENTS[component]
    frame = pd.read_csv(path, skiprows=find_header(path))
    frame.columns = [str(c).strip() for c in frame.columns]
    needed = ("Date", CUMULATIVE_COLUMN, UNCERTAINTY_COLUMN)
    if any(c not in frame.columns for c in needed):
        raise ValueError(f"{path.name}: missing required columns")
    frame = frame.loc[:, needed].copy()
    frame["Date"] = pd.to_datetime(frame["Date"])
    frame = frame.sort_values("Date").drop_duplicates("Date").reset_index(drop=True)
    frame = frame.loc[
        frame["Date"].dt.year.between(FIT_START_YEAR, LAST_OBS_YEAR)
    ].copy()
    frame[f"{component}_unrebased_cm"] = (
        -pd.to_numeric(frame[CUMULATIVE_COLUMN], errors="raise") * MM_TO_CM
    )
    frame[f"{component}_sigma_cm"] = (
        pd.to_numeric(frame[UNCERTAINTY_COLUMN], errors="raise") * MM_TO_CM
    )
    return frame[["Date", f"{component}_unrebased_cm", f"{component}_sigma_cm"]]


def calculate_constraints(source_dir):
    raw = read_component(source_dir, "E")
    for component in ("W", "G"):
        raw = raw.merge(
            read_component(source_dir, component),
            on="Date",
            how="inner",
            validate="one_to_one",
        )
    if (
        raw.empty
        or raw["Date"].iloc[0] != pd.Timestamp("1979-01-01")
        or raw["Date"].iloc[-1] != pd.Timestamp("2023-12-01")
    ):
        raise ValueError("Common E/W/G monthly record must span Jan 1979–Dec 2023")
    ref = raw["Date"].dt.year.between(REFERENCE_START_YEAR, REFERENCE_END_YEAR)
    for component in COMPONENTS:
        baseline = float(raw.loc[ref, f"{component}_unrebased_cm"].mean())
        raw[f"{component}_cm"] = raw[f"{component}_unrebased_cm"] - baseline
    raw["decimal_year"] = decimal_year(raw["Date"])
    raw["calendar_year"] = raw["Date"].dt.year

    # The original fitting uses calendar-year MEANS of monthly observations,
    # then averages COMPLETE four-year blocks beginning in 1979.
    annual_mean = raw.groupby("calendar_year", as_index=False).agg(
        decimal_year=("decimal_year", "mean"),
        **{f"{v}_cm": (f"{v}_cm", "mean") for v in COMPONENTS},
    )
    blocks = annual_mean.loc[
        annual_mean["calendar_year"].between(FIT_START_YEAR, FIT_END_YEAR)
    ].copy()
    blocks["block"] = (blocks["calendar_year"] - FIT_START_YEAR) // 4
    counts = blocks.groupby("block")["calendar_year"].transform("count")
    blocks = blocks.loc[counts.eq(4)]
    four = blocks.groupby("block", as_index=False).agg(
        decimal_year=("decimal_year", "mean"),
        **{f"{v}_cm": (f"{v}_cm", "mean") for v in COMPONENTS},
    )
    x2023 = float(raw["decimal_year"].iloc[-1])
    dt = four["decimal_year"].to_numpy(float) - x2023
    X = np.column_stack([np.ones_like(dt), dt, 0.5 * dt**2])
    Y = four[[f"{v}_cm" for v in COMPONENTS]].to_numpy(float)
    beta = np.linalg.lstsq(X, Y, rcond=None)[0]
    decembers = (
        raw.groupby("calendar_year", as_index=False)
        .tail(1)
        .sort_values("calendar_year")
        .reset_index(drop=True)
    )
    last24 = decembers.loc[decembers["calendar_year"].between(2000, 2023)].copy()
    years = last24["calendar_year"].to_numpy(dtype=np.int32)
    if not np.array_equal(years, np.arange(2000, 2024)):
        raise ValueError("Expected December data for each year 2000–2023")

    output = {"years": years, "boundary_decimal_year": np.asarray(x2023)}
    for idx, v in enumerate(COMPONENTS):
        observed_rate = float(
            np.gradient(
                decembers[f"{v}_cm"].to_numpy(float),
                decembers["decimal_year"].to_numpy(float),
                edge_order=2,
            )[-1]
        )
        yr2005 = raw["Date"].dt.year.eq(2005)
        if int(yr2005.sum()) != 12:
            raise ValueError("Expected 12 monthly observations in 2005")
        output[f"{v}_historical_cm"] = last24[f"{v}_cm"].to_numpy(dtype=np.float64)
        output[f"{v}_boundary_cm"] = np.asarray(float(beta[0, idx]))
        output[f"{v}_boundary_rate_cm_yr"] = np.asarray(observed_rate)
        output[f"{v}_reference_shift_cm"] = np.asarray(
            float(raw.loc[yr2005, f"{v}_cm"].mean())
        )
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-dir",
        type=Path,
        required=True,
        help="Folder containing the three original IMBIE CSVs",
    )
    parser.add_argument("--output", type=Path, required=True, help="Output .npz file")
    args = parser.parse_args()
    payload = calculate_constraints(args.input_dir)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.output, **payload)
    meta = {
        "description": "Deterministic IMBIE v3 historical constraints for SEJ2025",
        "reference_period": "1995-2014",
        "fitting_period": "1979-2022",
        "historical_output_years": "2000-2023 December observations",
        "boundary_source": "Quadratic fit to complete four-year means; observed December endpoint rate",
        "unit": "cm sea-level contribution",
        "source_sha256": {
            name: hashlib.sha256((args.input_dir / name).read_bytes()).hexdigest()
            for name in COMPONENTS.values()
        },
        "output_file": args.output.name,
        "output_sha256": hashlib.sha256(args.output.read_bytes()).hexdigest(),
    }
    meta_path = args.output.with_suffix(".json")
    meta_path.write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {args.output} and {meta_path}")


if __name__ == "__main__":
    main()
