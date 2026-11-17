"""Read immutable, precomputed IMBIE historical constraints for SEJ2025.

The corresponding offline generator is scripts/preprocessing/preprocess_imbie.py.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd

COMPONENTS = ('E', 'W', 'G')
EXPECTED_YEARS = np.arange(2000, 2024, dtype=np.int32)


def load_historical_constraints(path: str | Path):
    """Return (historical December values DataFrame, boundary dict, shifts dict).

    Contract matches the previously validated `interpolate_annual` inputs.
    Do not recalculate historical fitting during a FACTS2 run.
    """
    with np.load(Path(path), allow_pickle=False) as data:
        required = {'years', 'boundary_decimal_year'}
        required.update(f'{v}_{suffix}' for v in COMPONENTS for suffix in
                        ('historical_cm', 'boundary_cm', 'boundary_rate_cm_yr', 'reference_shift_cm'))
        missing = required - set(data.files)
        if missing:
            raise ValueError(f'Historical constraints missing: {sorted(missing)}')
        years = np.asarray(data['years'], dtype=int)
        if not np.array_equal(years, EXPECTED_YEARS):
            raise ValueError('Expected December observations for every year 2000–2023')
        boundary_date = float(data['boundary_decimal_year'])
        if not (2023 < boundary_date < 2024):
            raise ValueError('Historical boundary must be December 2023')
        historical = pd.DataFrame({'calendar_year': years})
        boundary = {}
        shifts = {}
        for v in COMPONENTS:
            history = np.asarray(data[f'{v}_historical_cm'], float)
            value = float(data[f'{v}_boundary_cm'])
            rate = float(data[f'{v}_boundary_rate_cm_yr'])
            shift = float(data[f'{v}_reference_shift_cm'])
            if history.shape != years.shape or not np.isfinite(history).all():
                raise ValueError(f'Invalid historical observations for {v}')
            if not np.isfinite([value, rate, shift]).all():
                raise ValueError(f'Non-finite boundary/reference values for {v}')
            historical[f'{v}_cm'] = history
            boundary[v] = (boundary_date, value, rate)
            shifts[v] = shift
    return historical, boundary, shifts
