"""Draw coupled 2050/2100/2200/2300 ice-sheet endpoints from a fixed CSV."""
from __future__ import annotations
import numpy as np
import pandas as pd
from .config import SCENARIOS, SCENARIO_CODES, CMIP6_SAMPLE_SEED, CMIP7_CRN_SEED


def crn_indices_by_sample_id(selection, n_rows, seed=CMIP7_CRN_SEED):
    """Shared FaIR IDs use the same row in a given scenario across pathways.

    """
    ids=pd.to_numeric(selection['sample_id'],errors='raise').to_numpy()
    if not np.isfinite(ids).all() or not np.array_equal(ids,np.floor(ids)):
        raise ValueError('CMIP7 CRN requires integer FaIR IDs')
    ids=ids.astype(np.int64)
    pathways=selection['ssp'].astype(str).to_numpy()
    sets=[]
    for p in pd.unique(pathways):
        group=ids[pathways==p]
        if len(np.unique(group))!=len(group): raise ValueError(f'Duplicate ID in {p}')
        sets.append(set(group.tolist()))
    if any(x!=sets[0] for x in sets[1:]):
        raise ValueError('All CMIP7 pathways must use identical FaIR sample-ID sets')
    ordered=np.array(sorted(sets[0]),dtype=np.int64)
    positions=np.searchsorted(ordered,ids)
    rng=np.random.default_rng(seed)
    streams={s:rng.integers(0,n_rows,size=len(ordered)) for s in SCENARIOS}
    selection_s=selection['selected_scenario'].to_numpy()
    result=np.full(len(selection),-1,dtype=np.int64)
    for s in SCENARIOS:
        mask=selection_s==s
        result[mask]=streams[s][positions[mask]]
    if (result<0).any():raise ValueError('Invalid selected warming scenario')
    return result


def sample_coupled_rows(selection, coupled, family='CMIP6', seed=None):
    """Select *whole rows*: never independently sample components or anchor years.

    Preserves CMIP6 with-replacement and CMIP7 sample-ID CRN behavior.
    """
    family=family.upper()
    if family=='CMIP6':
        rng=np.random.default_rng(CMIP6_SAMPLE_SEED if seed is None else seed)
        row_idx=rng.choice(len(coupled),size=len(selection),replace=True)
    elif family=='CMIP7':
        row_idx=crn_indices_by_sample_id(selection,len(coupled),CMIP7_CRN_SEED if seed is None else seed)
    else:raise ValueError('family must be CMIP6 or CMIP7')
    result=selection.reset_index(drop=True).copy()
    result['sej2025_row_index']=row_idx
    result['sej_sampling_rule']=('CMIP6_independent_rng_with_replacement' if family=='CMIP6'
                                 else 'CMIP7_CRN_shared_FaIR_sample_id_separate_warming_streams')
    result['coupled_trajectory_id']=(coupled.iloc[row_idx]['coupled_trajectory_id'].to_numpy()
                                      if 'coupled_trajectory_id' in coupled else row_idx)
    selected=result['selected_scenario'].to_numpy()
    for component in ('E','W','G'):
        result[f'KNOT_{component}_2050_cm']=np.nan
        for yr in ('21','22','23'):
            result[f'{component}_{yr}']=np.nan
    for scenario,code in SCENARIO_CODES.items():
        mask=selected==scenario
        if not mask.any():continue
        rows=coupled.iloc[row_idx[mask]]
        for component in ('E','W','G'):
            source=f'EMU2050_{component}{code}_cm'
            if source not in coupled:raise KeyError(source)
            result.loc[mask,f'KNOT_{component}_2050_cm']=rows[source].to_numpy(float)
            for year in ('21','22','23'):
                source=f'{component}{code}_{year}'
                if source not in coupled:raise KeyError(source)
                result.loc[mask,f'{component}_{year}']=rows[source].to_numpy(float)
    endpoints=[f'KNOT_{v}_2050_cm' for v in ('E','W','G')]+[f'{v}_{y}' for v in ('E','W','G') for y in ('21','22','23')]
    if result[endpoints].isna().any().any():raise ValueError('Missing sampled endpoint')
    return result
