import time

import numpy as np
import xarray as xr


def make_gmsl_projection_ds(
    global_samps: np.ndarray,
    years: np.ndarray,
    scenario: str,
    baseyear: int,
    module_name: str,
):
    """
    Create an xarray Dataset for global sea level rise projections from ice sheet samples.

    Parameters
    ----------
    global_samps : array-like
        Array of global sea level change samples, shape (samples, years).
    years : array-like
        Array of projection years.
    samples : array-like
        Array of sample indices.
    scenario : str
        Emissions scenario name.
    baseyear : int
        Reference year for projections.
    module_name : str
        Name of the module generating the projections.

    Returns
    -------
    xarray.Dataset
        Dataset containing sea level change data and metadata, with dimensions ('samples', 'years', 'locations').

    Notes
    -----
    - The returned dataset includes variables for sea level change, latitude, and longitude, as well as metadata attributes.
    - Latitude and longitude are set to NaN (np.inf) for global projections.
    """
    data = np.asarray(global_samps, dtype=np.float32)[:, :]
    nsamps = data.shape[0]
    ds = xr.Dataset(
        data_vars={
            "sea_level_change": (
                ("samples", "year"),
                data,
                {"units": "mm"},
            ),
        },
        coords={
            "year": (("year",), years),
            "samples": (("samples",), np.arange(nsamps)),
        },
        attrs={
            "description": f"Global SLR contribution from {module_name} module",
            "history": "Created " + time.ctime(time.time()),
            "scenario": scenario,
            "baseyear": baseyear,
        },
    )
    return ds
