"""Standalone integration runner for upstream GSAT, coupled SEJ, and historical IMBIE.

NOTE!!! The final FACTS2 module wrapper must still be matched to the framework's
experiment configuration and spatial fingerprint interface.
"""
import argparse
from pathlib import Path
import pandas as pd
from .scenario_selection import read_fair_gsat, select_scenarios
from .sampling import sample_coupled_rows
from .interpolation import interpolate_annual
from .historical import load_historical_constraints
from .config import YEARS
from .io import write_annual_ice_sheet_projections, write_component_gmsl_projections


# Build out your modules CLI entrypoint here
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--climate-data-file', '--gsat', dest='gsat', type=Path, required=True)
    parser.add_argument('--coupled',type=Path,required=True)
    parser.add_argument('--historical-constraints',type=Path,required=True)
    parser.add_argument('--family',choices=('CMIP6','CMIP7'),default='CMIP6')
    parser.add_argument('--scenario', '--pathway', dest='pathway', required=True)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--baseyear', type=int, default=2005,
                        help='Reference year for component GMSL outputs (default: 2005).')
    for component in ('EAIS', 'WAIS', 'GIS', 'AIS'):
        parser.add_argument(
            f'--output-{component}-gslr-file',
            type=Path,
            help=(
                f'Output path for {component} global NetCDF projections. '
                + (
                    'Aggregate diagnostic; do not pass to grd-fingerprinting.'
                    if component == 'AIS'
                    else 'Compatible with grd-fingerprinting.'
                )
            ),
        )
    args=parser.parse_args()
    temperature,ids=read_fair_gsat(args.gsat,family=args.family)
    selection=select_scenarios(temperature,ids,pathway=args.pathway)
    if args.family=='CMIP7':
        parser.error('CMIP7 common-random-number pairing requires all pathways in one batch; '
                     'use select_pathways() and sample_coupled_rows() together.')
    coupled=pd.read_csv(args.coupled)
    sampled=sample_coupled_rows(selection,coupled,family=args.family)
    hist,boundary,shifts=load_historical_constraints(args.historical_constraints)
    arrays,validation=interpolate_annual(sampled,hist,boundary,shifts)
    args.out.mkdir(parents=True,exist_ok=True)
    sampled.to_csv(args.out/'sampled_sej2025_realizations.csv',index=False)
    validation.to_csv(args.out/'trajectory_validation.csv',index=False)
    write_annual_ice_sheet_projections(
        arrays=arrays,
        sample_ids=ids,
        years=YEARS,
        output_file=args.out / 'annual_ice_sheet_trajectories_2000_2300.nc',
    )
    gmsl_outputs = {
        component: getattr(args, f'output_{component}_gslr_file')
        for component in ('EAIS', 'WAIS', 'GIS', 'AIS')
    }
    if any(path is not None for path in gmsl_outputs.values()):
        write_component_gmsl_projections(
            arrays=arrays,
            years=YEARS,
            sample_ids=ids,
            scenario=args.pathway,
            baseyear=args.baseyear,
            output_files=gmsl_outputs,
        )
    print('Saved:',args.out)

if __name__=='__main__':main()
