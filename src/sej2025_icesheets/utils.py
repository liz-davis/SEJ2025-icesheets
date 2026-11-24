"""Dataset construction helpers for SEJ2025 FACTS2 outputs."""

from __future__ import annotations

import time
import numpy as np
import xarray as xr


def make_gmsl_projection_ds(
    global_samps: np.ndarray,
    years: np.ndarray,
    scenario: str,
    baseyear: int,
    module_name: str,
    *,
    sample_ids=None,
    ice_source: str | None = None,
    input_reference_period: str = "1995-2014",
) -> xr.Dataset:
    """Create a global projection dataset.

    Parameters
    ----------
    global_samps : array-like
        Sea-level contribution in **centimeters** relative to
        ``input_reference_period``, shape (samples, years).
    years : array-like
        Projection years (must include ``baseyear``).
    scenario : str
        Emissions pathway label, e.g. 'ssp245'.
    baseyear : int
        Reference year for the *output* anomalies.
    module_name : str
        Name of generating module.
    sample_ids : array-like, optional
        Upstream FaIR sample identifiers in the supplied array order. Defaults
        to 0..N-1 for backward compatibility.
    ice_source : str, optional
        Component name for output metadata (EAIS, WAIS, GIS, or AIS).
    input_reference_period : str
        Baseline of the original series. Rebased outputs are relative to
        ``baseyear``; the original baseline is recorded as provenance.

    Returns
    -------
    xr.Dataset
        One variable ``sea_level_change`` in mm with dimensions
        (sample, year, budget_term, source_region).
    """
    values = np.asarray(global_samps)
    years = np.asarray(years)
    if values.ndim != 2:
        raise ValueError(f"Expected (samples, years), got {values.shape}")
    if years.ndim != 1 or values.shape[1] != len(years):
        raise ValueError("Year coordinate must match the second data axis")
    if not np.issubdtype(years.dtype, np.integer):
        raise ValueError("Years must be integer calendar years")
    if len(np.unique(years)) != len(years):
        raise ValueError("Year coordinate must contain unique years")
    if not np.isfinite(values).all():
        raise ValueError("GMSL projection arrays must contain only finite values")
    base_index = np.flatnonzero(years == baseyear)
    if len(base_index) != 1:
        raise ValueError(f"Requested baseyear {baseyear} not found in projection years")

    nsamps = values.shape[0]
    if sample_ids is None:
        sample_ids = np.arange(nsamps, dtype=np.int64)
    sample_ids = np.asarray(sample_ids)
    if sample_ids.ndim != 1 or len(sample_ids) != nsamps:
        raise ValueError("Sample identifiers must match first data axis")
    if len(np.unique(sample_ids)) != nsamps:
        raise ValueError("Sample identifiers must be unique")

    component_metadata = {
        "EAIS": ("ais", "EAIS"),
        "WAIS": ("ais", "WAIS"),
        "GIS": ("gris", "ALL"),
        "AIS": ("his", "AIS"), # Diagnostic only; exclude from GRD to avoid double counting
    }

    if ice_source not in component_metadata:
        raise ValueError(f"Unknown ice-sheet source: {ice_source}")

    budget_term, source_region = component_metadata[ice_source]

    # Rebase to the requested year and convert cm to mm.
    # Preserve float64 precision for the GRD interface.
    values = values.astype(np.float64)
 
    change_mm = (
        values - values[:, base_index[0], None]
    ) * 10.0

    ds = xr.Dataset(
        data_vars={
            "sea_level_change": (
                ("sample", "year", "budget_term", "source_region"),
                change_mm[:, :, None, None],
                {"units": "mm"},
            ),
        },
        coords={
            "sample": sample_ids,
            "year": years.astype(np.int32),
            "budget_term": [budget_term],
            "source_region": [source_region],
        },
        attrs={
            "Conventions": "CF-1.8",
            "title": "SEJ2025 Ice-Sheet GMSL Projections",
            "source": module_name,
            "history": "Created " + time.ctime(),
            "scenario": scenario,
            "baseyear": int(baseyear),
            "input_reference_period": input_reference_period,
        },
    )

    return ds