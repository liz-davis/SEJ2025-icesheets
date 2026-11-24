"""
Read and write projections generated within FACTS experiments.
 
This module holds the input/output helpers used at the different steps of a
FACTS experiment: the climate step, the global mean sea level (GMSL) step, and
the relative sea level (RSL) step. Each step has a paired reader and writer so
that on-disk conventions (NetCDF4 groups, variable names, compression settings)
live in one place.
 
Notes
-----
All writers serialise to NetCDF4. The climate-step writer additionally
validates the structure of the dataset and applies per-variable compression;
the GMSL and RSL writers pass the dataset through unmodified.
"""
import numpy as np
import xarray as xr
from typing import Literal

## The following functions handle IO for projections generated at different steps of a FACTS experiment

# Read & write climate step projections
def write_climate_step_projections(
        projections: xr.Dataset,
        output_file: str,
        output_type: Literal["gsat","ohc","deep_ocean_temperature", "climate"]
):
    """This is a function provided to write outputs from a FACTS climate module to a NetCDF4 file.
    Provide the prepared xr.Datset object that holds the projections, 
    the output path where the file will be written, 
    and specify the type of output """
    # Check expected dimensions, coords and data variables
    assert set(projections.dims) == {"years","samples"}, (
        f"Expected dimensions 'years' and 'samples', got {set(projections.dims)}"
    )
    assert set(projections.coords) == {"years","samples"},(
        f"Expected coordinates 'years' and 'samples', got {set(projections.coords)}"
    )

    expected_vars = {
        "gsat": {"surface_temperature"},
        "ohc": {"ohc"},
        "deep_ocean_temperature":{"deep_ocean_temperature"},
        "climate":{"surface_temperature","ocean_heat_content","deep_ocean_temperature"},
    }
    if output_type not in expected_vars:
        raise ValueError(
            f"Invalid output_type '{output_type}', expected one of 'gsat','ohc','deep_ocean_temperature', or 'climate'"
        )
    assert set(projections.data_vars) == expected_vars[output_type], (
        f"Expected data variable(s) {expected_vars[output_type]}, got {set(projections.data_vars)}"
    )
    encoding_by_type = {
        "gsat":{"surface_temperature": {"dtype": "float32","zlib": True, "complevel": 4}},
        "ohc":{"ohc": {"dtype":"float32","zlib":True,"complevel":4}},
        "deep_ocean_temperature":{
            "deep_ocean_temperature":{
            "dtype":"float32","zlib":True,"complevel":4}
            },
        "climate":{
            "surface_temperature":{"dtype":"float32","zlib":True, "complevel":4},
            "ocean_heat_content":{"dtype":"float32","zlib":True, "complevel":4},
            "deep_ocean_temperature":{"dtype":"float32","zlib":True, "complevel":4}
        }
    }
    return projections.to_netcdf(
        output_file,
        encoding = encoding_by_type[output_type]
    )

def read_climate(path: str):
    """
    Read the output of the climate step into memory.
 
    Parameters
    ----------
    path : str
        Path to the climate data file written by the FACTS climate step.

    Returns
    -------
    xarray.Dataset
        Dataset holding the climate projections.
 
    See Also
    --------
    write_climate_step_projections : Write climate-step projections to disk.
 
    Notes
    -----
    The dataset is opened lazily, so variable values are read from disk only
    when they are accessed. The underlying file handle stays open until the
    dataset is closed.
 
    Examples
    --------
    >>> ds = read_climate("climate.nc")  # doctest: +SKIP
    """
    return xr.open_dataset(path)


# Read and write gmsl projections
def write_gmsl_projections(
    projections: xr.Dataset,
    output_file: str,
):
    """Write GRD-compatible GMSL projections to NetCDF4."""

    encoding = {
        "sea_level_change": {
            "dtype": "float64",
            "zlib": True,
            "complevel": 4,
        }
    }

    projections.to_netcdf(
        output_file,
        format="NETCDF4",
        engine="netcdf4",
        encoding=encoding,
    )


def read_gmsl_projections(path: str)-> xr.Dataset:
    """
    Read global mean sea level projections into memory.
 
    Parameters
    ----------
    path : str
        Path to the NetCDF file holding the GMSL projections.
 
    Returns
    -------
    xarray.Dataset
        Dataset holding the GMSL projections.
 
    See Also
    --------
    write_gmsl_projections : Write GMSL projections to disk.
 
    Notes
    -----
    The dataset is opened lazily, so variable values are read from disk only
    when they are accessed.
    """
    return xr.open_dataset(path)

# Read and write rsl projections
def write_rsl_projections(projections:xr.Dataset, output_file:str):
    """
    Write relative sea level projections to a NetCDF file.
 
    Parameters
    ----------
    projections : xarray.Dataset
        Dataset holding the RSL projections to serialise.
    output_file : str
        Path of the NetCDF file to write. An existing file is overwritten.
 
    Returns
    -------
    None

    See Also
    --------
    read_rsl_projections : Read RSL projections back into memory.
 
    Notes
    -----
    The dataset is written as-is: no validation is performed and no encoding
    or compression is applied.
    """
    projections.to_netcdf(output_file)

def read_rsl_projections(path: str)->xr.Dataset:
    """
    Read relative sea level projections into memory.
 
    Parameters
    ----------
    path : str
        Path to the NetCDF file holding the RSL projections.
 
    Returns
    -------
    xarray.Dataset
        Dataset holding the RSL projections.
 
    See Also
    --------
    write_rsl_projections : Write RSL projections to disk.
 
    Notes
    -----
    The dataset is opened lazily, so variable values are read from disk only
    when they are accessed.
    """
    return xr.open_dataset(path)

#---- Add your own functions to read/write additional data here ----#


# SEJ2025 annual ice-sheet projection output
def write_annual_ice_sheet_projections(
    arrays: dict,
    sample_ids,
    years,
    output_file: str,
):
    """Write annual EAIS, WAIS, GrIS, AIS, and total ice-sheet contributions.

    Arrays must have shape (samples, years), in centimeters relative to
    the 1995–2014 reference period.
    """
    sample_ids = np.asarray(sample_ids)
    years = np.asarray(years)
    if sample_ids.ndim != 1 or years.ndim != 1:
        raise ValueError("Sample IDs and years must be one-dimensional.")
    if len(np.unique(sample_ids)) != len(sample_ids):
        raise ValueError("Sample identifiers must be unique.")
    if len(np.unique(years)) != len(years):
        raise ValueError("Years must be unique.")
    if not arrays:
        raise ValueError("No annual projection arrays provided.")

    expected_shape = (len(sample_ids), len(years))
    dataset_vars = {}
    for name, values in arrays.items():
        values = np.asarray(values)
        if values.shape != expected_shape:
            raise ValueError(
                f"{name} has shape {values.shape}; expected {expected_shape}."
            )
        dataset_vars[name] = (("samples", "years"), values)

    ds = xr.Dataset(
        dataset_vars,
        coords={"samples": sample_ids, "years": years},
    )
    ds.attrs.update({
        "title": "SEJ2025 Annual Ice-Sheet Projections",
        "description": (
            "Annual ice-sheet contributions from temporally coupled SEJ2025 "
            "Monte Carlo realizations, emulandice 2050, and historical IMBIE."
        ),
        "reference_period": "1995-2014",
    })
    for name in arrays:
        ds[name].attrs["units"] = "cm"

    encoding = {
        name: {"dtype": "float32", "zlib": True, "complevel": 4}
        for name in arrays
    }
    return ds.to_netcdf(
        output_file, format="NETCDF4", engine="netcdf4", encoding=encoding
    )


def write_component_gmsl_projections(
    arrays: dict,
    years,
    sample_ids,
    scenario: str,
    baseyear: int,
    output_files: dict,
):
    """Write GRD-compatible EAIS, WAIS, and GIS GMSL projections. 
       AIS is supported as an optional diagnostic output only.

    ``arrays`` are SEJ annual trajectories in cm, keys E/W/G/AIS.
    ``output_files`` uses component keys EAIS/WAIS/GIS/AIS; a None path skips
    that component. Outputs are mm relative to ``baseyear``. This does not
    modify the input arrays or the combined annual scientific file.
    """
    from pathlib import Path
    from .utils import make_gmsl_projection_ds

    source_key = {"EAIS": "E", "WAIS": "W", "GIS": "G", "AIS": "AIS"}
    unknown = set(output_files) - set(source_key)
    if unknown:
        raise ValueError(f"Unknown GMSL output component(s): {sorted(unknown)}")
    for component, path in output_files.items():
        if path is None:
            continue
        if source_key[component] not in arrays:
            raise ValueError(f"Missing annual component {source_key[component]} for {component}")
        ds = make_gmsl_projection_ds(
            global_samps=arrays[source_key[component]],
            years=years,
            sample_ids=sample_ids,
            scenario=scenario,
            baseyear=baseyear,
            module_name="SEJ2025-icesheets",
            ice_source=component,
        )
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        write_gmsl_projections(ds, str(target))
