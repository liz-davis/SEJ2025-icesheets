# SEJ25 emulandice–SEJ temporally coupled ensemble

**Data file:** `SEJ25_emulandice2050_recoupled_by_warming.csv`  
**Purpose:** Precomputed joint realizations for the 2025 structured-expert-judgement (SEJ) ice-sheet projections, with emulandice-derived 2050 constraints and SEJ projections at 2100, 2200, and 2300. This file is an input to subsequent climate-scenario selection, sampling, and annual trajectory generation; **it is not itself an annual projection file**.

## File summary

| Property | Value |
| --- | --- |
| File format | CSV (header row; one realization per row) |
| Realizations | 50,000 |
| Columns | 208 |
| Realization identifier | `coupled_trajectory_id` (unique integers 0–49,999) |
| Warming scenarios | LOW (`2`, 2°C), INTERMEDIATE (`35`, 3.5°C), HIGH (`5`, 5°C) |
| Projection anchor years | 2050, 2100, 2200, 2300 |
| Ice-sheet components | East Antarctic Ice Sheet (EAIS), West Antarctic Ice Sheet (WAIS), Greenland Ice Sheet (GrIS) |

## Provenance and construction

The dataset combines (1) 2050 ice-sheet projections informed by **emulandice** and (2) the 2025 SEJ distributions at 2100, 2200, and 2300. These inputs were coupled offline to impose the selected within-scenario temporal dependence structure across anchor years and ice-sheet components. The temporal-coupling preprocessing, including the treatment of the 2050 emulandice anchor, is maintained separately from the downstream projection generation code.


Relevant upstream script:

- [`02_TEMPORAL_COUPLING_EMULANDICE2050_12VAR_MAXDET.py`](https://github.com/bobkopp/ice-sheet-sej-2025/blob/main/scripts/interpolation/temporal_dependence/02_TEMPORAL_COUPLING_EMULANDICE2050_12VAR_MAXDET.py)

**Important:** The three warming scenarios are separate scenario-conditional ensembles. A common row number across different warming scenarios does *not*, by itself, imply a physically matched realization across warming levels. The documented downstream workflow selects **one** scenario for a FaIR member and extracts all anchor-year/component values from **one row within that selected scenario**.

## Column naming conventions

### SEJ projections at 2100, 2200, and 2300

Columns follow the general pattern:

```text
<variable><warming_code>_<year_code>
```

- `warming_code`: `2` = LOW (2°C), `35` = INTERMEDIATE (3.5°C), `5` = HIGH (5°C).
- `year_code`: `21` = 2100, `22` = 2200, `23` = 2300.
- The component totals used by the annual trajectory engine are `E` (EAIS), `W` (WAIS), and `G` (GrIS).

For example, `E2_21`, `W2_21`, and `G2_21` are the EAIS, WAIS, and GrIS contributions in 2100 for the LOW warming scenario. Corresponding fields exist at 2200 and 2300 and for the other warming scenarios.

The file also includes more detailed SEJ component fields, for example `EA`, `ED`, `ER`, `WA`, `WD`, `WR`, `GA`, `GD`, and `GR`, as well as `SLR` aggregate fields. These additional fields are retained from the original coupled dataset; their detailed scientific definitions should be taken from the upstream SEJ data dictionary or processing scripts rather than inferred from the abbreviated column names alone.

**Units:** The downstream annual-projection script treats SEJ endpoint fields such as `E2_21`, `W2_21`, and `G2_21` as **millimeters**, multiplies them by `0.1` to convert to centimeters, and applies its reference-period adjustments. Keep these conversions when reusing the dataset.

### emulandice 2050 values

Columns include:

```text
EMU2050_E2_cm    EMU2050_W2_cm    EMU2050_G2_cm    EMU2050_AIS2_cm
EMU2050_E35_cm   EMU2050_W35_cm   EMU2050_G35_cm   EMU2050_AIS35_cm
EMU2050_E5_cm    EMU2050_W5_cm    EMU2050_G5_cm    EMU2050_AIS5_cm
```

All of these fields are already in **centimeters**. Each warming level also has an `EMU2050_PAIR_ID_<warming_code>` column to retain the emulandice pairing identifier.

The 2050 values are **already attached to the coupled realization**: no independent 2050 emulandice resampling or rank matching is required during downstream FACTS2 execution.

### Observation and diagnostic fields

Each warming level has fields such as:

- `OBS2_E_fitted_2023_cm`, `OBS2_W_fitted_2023_cm`, `OBS2_G_fitted_2023_cm`
- `OBS2_E_rate_cm_per_yr`, `OBS2_E_acceleration_cm_per_yr2`
- `OBS2_E_2050_cm`, `OBS2_E_rank_2050`

Equivalent `OBS35_*` and `OBS5_*` columns are provided. These preserve observation-based quantities and rank information used during preprocessing. They should not be confused with the **final** `EMU2050_*_cm` anchor values selected by the downstream annual trajectory engine.

## How the existing projection workflow uses this file

1. **Select a warming scenario for each FaIR realization.** The scenario-selection scripts compare cumulative FaIR warming at **2100 and 2200** to the LOW, INTERMEDIATE, and HIGH SEJ reference warming trajectories, using range-normalized Euclidean distance and inverse-distance weights (`p = 1`). They force LOW or HIGH only when the realization lies beyond the corresponding endpoint reference **at both checkpoints in the same direction**. One selected warming scenario applies to the full ice-sheet trajectory.
2. **Sample a complete coupled row.** `sample_coupled_rows()` in `core_annual_tension.py` chooses a row index within the selected warming-scenario ensemble and extracts the corresponding `E`, `W`, and `G` endpoints for 2100/2200/2300 **and** the `EMU2050_*_cm` anchor values from that same row.
3. **Build annual trajectories.** The annual engine uses IMBIE-based historical context through 2023 and constructs each component's 2024–2300 projection through the fitted 2023 boundary, the recoupled 2050 emulandice anchor, and the three later SEJ anchors using a classical spline under tension (`p = 0.03` in the supplied implementation).
4. **Derive total contributions.** Annual total Antarctic and total ice-sheet contributions are calculated from the component arrays, rather than requiring separate independently sampled total trajectories.


## Use in FACTS2

This CSV is intended to be a **read-only, versioned input** in the new `SEJ2025-icesheets` FACTS2 module (for example under `module_specific_input_data/`). A FACTS2 run should perform scenario selection, row sampling, annual interpolation, and output formatting, **not** rerun the offline emulandice/SEJ temporal-coupling optimization.


## Reproducibility and citation

Source code repository: https://github.com/bobkopp/ice-sheet-sej-2025

---

*This README documents the structure of the supplied CSV and the behavior visible in the accompanying scenario-selection and annual-projection scripts. It is not a substitute for the complete scientific methods or upstream column-level data dictionary.*
