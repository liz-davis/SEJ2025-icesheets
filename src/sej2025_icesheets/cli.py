"""
Standalone integration runner for upstream GSAT, coupled SEJ, and historical IMBIE.

NOTE!!! The final FACTS2 module wrapper must still be matched to the framework's
experiment configuration and spatial fingerprint interface.
"""

import click
import argparse
from pathlib import Path
import pandas as pd
import xarray as xr
from .scenario_selection import read_fair_gsat, select_scenarios
from .sampling import sample_coupled_rows
from .interpolation import historical_context, interpolate_annual
from .config import YEARS

# Build out your modules CLI entrypoint here
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gsat',type=Path,required=True)
    parser.add_argument('--coupled',type=Path,required=True)
    parser.add_argument('--historical-project-dir',type=Path,required=True)
    parser.add_argument('--family',choices=('CMIP6','CMIP7'),default='CMIP6')
    parser.add_argument('--pathway',required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    temperature,ids=read_fair_gsat(args.gsat,family=args.family)
    selection=select_scenarios(temperature,ids,pathway=args.pathway)
    if args.family=='CMIP7':
        parser.error('CMIP7 common-random-number pairing requires all pathways in one batch; '
                     'use select_pathways() and sample_coupled_rows() together.')
    coupled=pd.read_csv(args.coupled)
    sampled=sample_coupled_rows(selection,coupled,family=args.family)
    hist,boundary,shifts=historical_context(args.historical_project_dir)
    arrays,validation=interpolate_annual(sampled,hist,boundary,shifts)
    args.out.mkdir(parents=True,exist_ok=True)
    sampled.to_csv(args.out/'sampled_sej2025_realizations.csv',index=False)
    validation.to_csv(args.out/'trajectory_validation.csv',index=False)
    ds=xr.Dataset({k:(('samples','years'),v) for k,v in arrays.items()},
                  coords={'samples':ids,'years':YEARS})
    ds.attrs['reference_period']='1995-2014'
    for k in arrays:ds[k].attrs['units']='cm'
    ds.to_netcdf(args.out/'annual_ice_sheet_trajectories_2000_2300.nc')
    print('Saved:',args.out)

if __name__=='__main__':main()
