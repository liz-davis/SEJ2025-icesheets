"""FaIR-to-SEJ scenario selection shared by CMIP6 and CMIP7.

Preserves the adopted cumulative-warming, two-checkpoint, inverse-distance
sampling method and endpoint truncation. Inputs are warming anomalies in K
on the same baseline as the reference warming trajectories.
"""
from __future__ import annotations
import json
import warnings
import numpy as np
import pandas as pd
import xarray as xr
from .config import (YEARS, SCENARIOS, REFERENCE_YEARS, REFERENCE_WARMING,
                     CHECKPOINTS, SELECTION_SEED)

METHOD = 'range_normalized_checkpoint_euclidean_idw_p1_double_checkpoint_truncation'
TRUNCATION_RULE = 'outside_same_direction_at_both_2100_and_2200'


def read_fair_gsat(path, family='CMIP6', variable=None, warming_offset=0.0):
    """Read upstream FACTS2 FaIR GSAT NetCDF and return (samples, years).

    Accepts (samples, years, locations) and selects global location -1.
    Requires an annual 2000-2300 subset and SEJ-compatible anomaly baseline.
    """
    with xr.open_dataset(path) as ds:
        if variable is None:
            candidates = [x for x in ('surface_temperature','gsat','GSAT',
                          'global_mean_surface_temperature') if x in ds]
            if len(candidates) != 1:
                raise ValueError(f'Cannot uniquely identify GSAT variable: {list(ds.data_vars)}')
            variable = candidates[0]
        centered = ds.attrs.get('Centered', ds.attrs.get('centered'))
        if isinstance(centered, bytes):
            centered = centered.decode('utf-8')
        if centered is not None:
            baseline = str(centered).lower().replace('–', '-').replace(' ', '')
            if '1850-1900' not in baseline:
                raise ValueError(
                    f'FaIR warming baseline {centered!r} does not match the '
                    'SEJ reference baseline (1850-1900).'
                )
        else:
            warnings.warn(
                'FaIR NetCDF has no Centered metadata; independently verify '
                'that GSAT is relative to 1850-1900.', stacklevel=2
            )
        da = ds[variable]
        year_dim = next((d for d in ('years','year','time') if d in da.dims), None)
        sample_dim = next((d for d in ('samples','sample','realization','realizations',
                                     'draw','ensemble_member') if d in da.dims), None)
        if year_dim is None or sample_dim is None:
            raise ValueError(f'Expected year/sample dimensions, got {da.dims}')
        for dim in tuple(da.dims):
            if dim in (year_dim, sample_dim):
                continue
            if dim == 'locations' and -1 in da[dim].values:
                da = da.sel({dim:-1})
            elif da.sizes[dim] == 1:
                da = da.isel({dim:0})
            else:
                raise ValueError(f'Ambiguous extra dimension {dim}')
        coord = np.asarray(da[year_dim].values)
        years = (pd.DatetimeIndex(coord).year.to_numpy() if np.issubdtype(coord.dtype, np.datetime64)
                 else coord.astype(int))
        if len(set(years)) != len(years):
            raise ValueError('Duplicate years')
        chosen = np.flatnonzero((years >= YEARS[0]) & (years <= YEARS[-1]))
        if not np.array_equal(years[chosen], YEARS):
            raise ValueError('FaIR must have complete annual 2000–2300 values')
        da = da.isel({year_dim:chosen}).transpose(sample_dim, year_dim)
        temperatures = np.asarray(da.values, dtype=float) + warming_offset
        sample_ids = np.asarray(da[sample_dim].values)
        if not np.isfinite(temperatures).all() or len(set(map(str,sample_ids))) != len(sample_ids):
            raise ValueError('Nonfinite GSAT or duplicate sample identifiers')
        median_2000 = np.median(temperatures[:, 0]); median_2020 = np.median(temperatures[:,20])
        if family.upper() == 'CMIP7':
            if not (-1 <= median_2000 <= 3 and -1 <= median_2020 <= 6):
                raise ValueError('CMIP7 GSAT appears to use an incompatible baseline')
            if abs(median_2000-.7)>.45 or abs(median_2020-1.2)>.65:
                warnings.warn('Check CMIP7 GSAT anomaly baseline against SEJ references', stacklevel=2)
    return temperatures, sample_ids


def cumulative_warming(temperature):
    temperature=np.asarray(temperature, dtype=float)
    if temperature.ndim != 2 or temperature.shape[1] != len(YEARS):
        raise ValueError('Expected [samples, 301 annual years]')
    out=np.zeros_like(temperature)
    out[:,1:]=np.cumsum(.5*(temperature[:,:-1]+temperature[:,1:]),axis=1)
    return out


def reference_cumulative():
    reference=np.stack([np.interp(YEARS,REFERENCE_YEARS,REFERENCE_WARMING[s]) for s in SCENARIOS])
    return cumulative_warming(reference)


def inverse_distance_weights(distances):
    distances=np.asarray(distances,dtype=float)
    if distances.ndim!=2 or distances.shape[1]!=3 or not np.isfinite(distances).all() or (distances<0).any():
        raise ValueError('distances must be finite nonnegative N x 3')
    weights=np.zeros_like(distances)
    matches=np.isclose(distances,0,rtol=0,atol=1e-12)
    exact=matches.any(axis=1)
    weights[exact]=matches[exact]/matches[exact].sum(axis=1,keepdims=True)
    inv=1/distances[~exact]
    weights[~exact]=inv/inv.sum(axis=1,keepdims=True)
    return weights


def select_scenarios(temperature, sample_ids=None, pathway='', rng=None):
    """Select ONE warming class for each FaIR member, retaining diagnostics.

    Supply a shared RNG to reproduce the original sequential multi-pathway
    draw stream; standalone runs default to the production selection seed.
    """
    if rng is None: rng=np.random.default_rng(SELECTION_SEED)
    c=cumulative_warming(temperature)
    reference=reference_cumulative()
    ki=np.array([np.where(YEARS==y)[0][0] for y in CHECKPOINTS])
    span=reference[2,ki]-reference[0,ki]
    distances=np.sqrt(np.sum(((c[:,None,ki]-reference[None,:,ki])/span[None,None,:])**2,axis=2))
    raw=inverse_distance_weights(distances)
    low=(c[:,ki]<=reference[0,ki]).all(axis=1)
    high=(c[:,ki]>=reference[2,ki]).all(axis=1)
    weights=raw.copy()
    weights[low]=[1.,0.,0.]; weights[high]=[0.,0.,1.]
    uniforms=rng.random(len(c))
    cdf=np.cumsum(weights,axis=1);cdf[:,-1]=1.
    selected=np.asarray(SCENARIOS,dtype=object)[(uniforms[:,None]>cdf).sum(axis=1)]
    if sample_ids is None: sample_ids=np.arange(len(c))
    result=pd.DataFrame({'ssp':pathway,'sample_id':sample_ids,'selected_scenario':selected,
        'sampling_method':METHOD,'selection_uniform':uniforms,'horizon':2300,
        'checkpoint_years':json.dumps(list(CHECKPOINTS)), 'truncation_years':json.dumps(list(CHECKPOINTS)),
        'truncation_rule':TRUNCATION_RULE,'idw_power':1.0,
        'selection_rule':np.where(low,'TRUNCATED_LOW_BOTH_CHECKPOINTS',
                                  np.where(high,'TRUNCATED_HIGH_BOTH_CHECKPOINTS','PROBABILISTIC_IDW')),
        'reference_range_class':np.where(low,'BELOW_LOW_BOTH',
                                         np.where(high,'ABOVE_HIGH_BOTH','NOT_DOUBLE_OUTSIDE'))})
    for i,y in enumerate(CHECKPOINTS):
        result[f'cumulative_warming_{y}_for_truncation']=c[:,ki[i]]
        result[f'reference_range_position_{y}']=(c[:,ki[i]]-reference[0,ki[i]])/span[i]
        result[f'below_low_{y}']=c[:,ki[i]]<=reference[0,ki[i]]
        result[f'above_high_{y}']=c[:,ki[i]]>=reference[2,ki[i]]
    result['reference_range_position']=.5*(result['reference_range_position_2100']+result['reference_range_position_2200'])
    for i,s in enumerate(('low','intermediate','high')):
        result[f'distance_{s}']=distances[:,i]
        result[f'raw_weight_{s}']=raw[:,i]
        result[f'weight_{s}']=weights[:,i]
    return result


def select_pathways(fair_paths, family='CMIP6', seed=SELECTION_SEED, **read_kwargs):
    """Process pathways in insertion order, sharing the original selection RNG."""
    rng=np.random.default_rng(seed)
    frames=[]
    for pathway,path in fair_paths.items():
        data,ids=read_fair_gsat(path,family=family,**read_kwargs)
        frames.append(select_scenarios(data,ids,pathway=pathway,rng=rng))
    return pd.concat(frames,ignore_index=True)
