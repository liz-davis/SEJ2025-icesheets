# Historical ice-sheet observations (IMBIE v3)

This directory contains historical ice-sheet mass-balance observations used by the **SEJ2025-icesheets** module of the Framework for Assessing Changes To Sea-level (FACTS2). The three comma-separated files come from the *Ice Sheet Mass Balance Intercomparison Exercise* (IMBIE) and provide historical estimates for **East Antarctica (EAIS)**, **West Antarctica (WAIS)**, and **Greenland (GrIS)** through December 2023.

These are **observational inputs, not SEJ Monte Carlo samples**. They anchor the historical portion of the ice-sheet sea-level contribution time series and provide the observed endpoint value and rate for the transition to the future trajectories constructed from the coupled 2050–2300 projections.

## Files

| Filename | Region | First monthly record | Last monthly record | Records |
| --- | --- | --- | --- | ---: |
| `imbie3_east_antarctica_mm_partitioned.csv` | East Antarctic Ice Sheet (EAIS) | 1979-01-01 | 2023-12-01 | 540 |
| `imbie3_west_antarctica_mm_partitioned.csv` | West Antarctic Ice Sheet (WAIS) | 1979-01-01 | 2023-12-01 | 540 |
| `imbie3_greenland_mm_partitioned.csv` | Greenland Ice Sheet (GrIS) | 1971-07-01 | 2023-12-01 | 630 |

The inventories above are based on the provided files. Each CSV contains **13 columns**: a date plus 12 quantities.

## Data provenance

The source dataset is **“Mass balance of the Greenland and Antarctic Ice Sheets from the 1970s to 2023,” version 1.0**, issued by the NERC Environmental Data Service / UK Polar Data Centre as part of IMBIE. Its reconciled estimates combine satellite altimetry, gravimetry, and input–output observations. The dataset metadata describe 45 independent estimates: 21 for Antarctica and 24 for Greenland.

- **Dataset DOI:** https://doi.org/10.5285/128c5e33-5224-4197-82f0-19dcc95b80a0
- **Data publication:** Otosaka et al. (2026), *Mass balance of the Greenland and Antarctic Ice Sheets from the 1970s to 2023*, version 1.0, NERC EDS UK Polar Data Centre.
- **Associated journal publication:** Otosaka et al. (2026), *Mass balance of the Greenland and Antarctic ice sheets from the 1970s to 2023*, *Scientific Data*, **13**, 1301. https://doi.org/10.1038/s41597-026-08088-0
- **IMBIE processing software:** Otosaka et al. (2024), *IMBIE software v3.1*, https://doi.org/10.5281/zenodo.14217193
- **Project:** https://imbie.org/
- **License recorded in the source file:** UK Open Government Licence, version 3.0 — https://www.nationalarchives.gov.uk/doc/open-government-licence/version/3/


## Column definitions

The same schema appears in all three CSVs:

| Column | Units | Interpretation |
| --- | --- | --- |
| `Date` | YYYY-MM-DD | Monthly observation date |
| `Mass balance (mm/yr)` | mm/yr | Total mass-balance rate expressed in sea-level-equivalent units |
| `Mass balance uncertainty (mm/yr)` | mm/yr | Reported uncertainty of the total rate |
| `Cumulative mass balance anomaly (mm)` | mm | Cumulative total mass-balance anomaly |
| `Cumulative mass balance anomaly uncertainty (mm)` | mm | Reported uncertainty of the cumulative total |
| `Surface mass balance anomaly (mm/yr)` | mm/yr | Surface-mass-balance contribution to the rate |
| `Surface mass balance anomaly uncertainty (mm/yr)` | mm/yr | Reported uncertainty of the surface component rate |
| `Cumulative surface mass balance anomaly (mm)` | mm | Cumulative surface component |
| `Cumulative surface mass balance anomaly uncertainty (mm)` | mm | Reported uncertainty of the cumulative surface component |
| `Dynamics mass balance anomaly (mm/yr)` | mm/yr | Ice-dynamics contribution to the rate |
| `Dynamics mass balance anomaly uncertainty (mm/yr)` | mm/yr | Reported uncertainty of the dynamics rate |
| `Cumulative dynamics mass balance anomaly (mm)` | mm | Cumulative dynamics component |
| `Cumulative dynamics mass balance anomaly uncertainty (mm)` | mm | Reported uncertainty of the cumulative dynamics component |

**Units and conventions:** The `mm` and `mm/yr` labels are sea-level-equivalent quantities. Preserve the source sign conventions; do not reverse signs or interpret a negative mass-balance anomaly as positive sea-level contribution without checking the relevant processing code. Uncertainty columns are described by the source metadata as **one-sigma** uncertainties. The source CSVs should remain unchanged; any unit conversions and reference-period adjustments belong in module code.

## Role in SEJ2025 annual projection construction

In the SEJ research workflow, the shared annual-trajectory implementation (`core_annual_tension.py`) uses historical IMBIE observations to construct the historical part of the projection and the boundary conditions for the future:

1. Read the partitioned IMBIE records for EAIS, WAIS, and GrIS using the historical-data helper functions.
2. Construct the historical output for **2000–2023**, rebased to the **1995–2014 reference period** used by that workflow.
3. Estimate the fitted end-of-2023 historical value and observed year-end rate for each component. These supply the **left endpoint value and slope** of the future spline.
4. Connect the historical boundary to four coupled future knots: **2050 (emulandice)** and **2100, 2200, 2300 (SEJ)**.
5. Generate annual values for **2024–2300** using the adopted classical spline under tension (parameter **p = 0.03**), maintaining continuity with the historical boundary.
6. Derive combined components (for example, total Antarctic and total ice-sheet contribution) from the component trajectories rather than treating them as independent historical inputs.

**Important:** The observational time series are not independent random historical draws for each SEJ realization. In the current research workflow, the same fitted historical boundary is used for the sampled future trajectories of each component. The input files also include published observational uncertainties, but these are not automatically converted into sampled uncertainty in that procedure.


## Suggested citation

> Otosaka, I. N., et al. (2026). *Mass balance of the Greenland and Antarctic Ice Sheets from the 1970s to 2023* (Version 1.0) [Data set]. NERC EDS UK Polar Data Centre. https://doi.org/10.5285/128c5e33-5224-4197-82f0-19dcc95b80a0

Please cite the dataset DOI when using these files, and the associated *Scientific Data* article when describing the reconciliation methodology.
