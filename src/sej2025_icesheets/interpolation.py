"""
Annual SEJ2025 sea-level trajectories using tension spline.

"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .config import TENSION_P, YEARS

MM_TO_CM = 0.1


class TensionSpline:
    """Classical interpolating spline under tension with continuous first derivative."""

    def __init__(self, x, y, p, left_slope):
        self.x = np.asarray(x, float)
        y = np.asarray(y, float)
        if y.ndim == 1:
            y = y[:, None]
            self.squeeze = True
        else:
            self.squeeze = False
        self.y = y
        self.p = float(p)
        slope = np.asarray(left_slope, float)
        if slope.ndim == 0:
            slope = np.full(y.shape[1], float(slope))
        self.left_slope = slope
        self.M = self._solve()

    def _solve(self):
        x, y, p = self.x, self.y, self.p
        h = np.diff(x)
        delta = np.diff(y, axis=0) / h[:, None]
        n = len(x)
        A = np.zeros((n - 1, n - 1), float)
        rhs = np.empty((n - 1, y.shape[1]), float)

        h0 = h[0]
        s0 = np.sinh(p * h0)
        c0 = np.cosh(p * h0) / s0
        A[0, 0] = c0 / p - 1.0 / (p**2 * h0)
        if n - 1 > 1:
            A[0, 1] = 1.0 / (p**2 * h0) - 1.0 / (p * s0)
        rhs[0] = delta[0] - self.left_slope

        for i in range(1, n - 1):
            hl, hr = h[i - 1], h[i]
            sl, sr = np.sinh(p * hl), np.sinh(p * hr)
            cl, cr = np.cosh(p * hl) / sl, np.cosh(p * hr) / sr
            A[i, i - 1] = 1.0 / (p**2 * hl) - 1.0 / (p * sl)
            A[i, i] = cl / p - 1.0 / (p**2 * hl) + cr / p - 1.0 / (p**2 * hr)
            if i + 1 < n - 1:
                A[i, i + 1] = 1.0 / (p**2 * hr) - 1.0 / (p * sr)
            rhs[i] = delta[i] - delta[i - 1]

        M = np.zeros_like(y)
        M[:-1] = np.linalg.solve(A, rhs)
        return M

    def _eval(self, xx, derivative=False):
        xx = np.atleast_1d(np.asarray(xx, float))
        out = np.empty((len(xx), self.y.shape[1]), float)
        idx = np.searchsorted(self.x, xx, side="right") - 1
        idx = np.clip(idx, 0, len(self.x) - 2)
        p = self.p
        for k in np.unique(idx):
            mask = idx == k
            u = xx[mask] - self.x[k]
            h = self.x[k + 1] - self.x[k]
            y0, y1 = self.y[k], self.y[k + 1]
            M0, M1 = self.M[k], self.M[k + 1]
            sh = np.sinh(p * h)
            if derivative:
                term0 = -np.cosh(p * (h - u))[:, None] / (p * sh)
                term1 = np.cosh(p * u)[:, None] / (p * sh)
                linear = (y1 - y0) / h + (M0 - M1) / (p**2 * h)
                out[mask] = term0 * M0 + term1 * M1 + linear
            else:
                lh = np.sinh(p * (h - u))[:, None] / (p**2 * sh)
                rh = np.sinh(p * u)[:, None] / (p**2 * sh)
                ll = ((h - u) / h)[:, None]
                rl = (u / h)[:, None]
                out[mask] = (
                    lh * M0 + rh * M1 + ll * (y0 - M0 / p**2) + rl * (y1 - M1 / p**2)
                )
        return out[:, 0] if self.squeeze else out

    def value(self, xx):
        return self._eval(xx, derivative=False)

    def rate(self, xx):
        return self._eval(xx, derivative=True)


# Historical fitting is performed offline by scripts/preprocessing/preprocess_imbie.py.
# Use `load_historical_constraints()` to supply immutable arrays to interpolation.


def interpolate_annual(sampled, annual_hist, boundary, shifts, tension_p=TENSION_P):
    """Output cm, shape (sample, year) for E, W, G, AIS, SLR.

    historical columns are E_cm, W_cm and G_cm, with calendar_year.
    The SEJ 2100/2200/2300 endpoint columns are in mm (2005 baseline);
    knots at 2050 and historical observations are in cm.
    """
    n = len(sampled)
    output = {}
    validations = []
    future = YEARS >= 2024
    for v in ("E", "W", "G"):
        full = np.full((n, len(YEARS)), np.nan, dtype=np.float32)
        for row in annual_hist.itertuples():
            yr = int(row.calendar_year)
            if 2000 <= yr <= 2023:
                full[:, yr - 2000] = np.float32(getattr(row, f"{v}_cm"))
        x2023, y2023, rate2023 = boundary[v]
        knots = {2050: sampled[f"KNOT_{v}_2050_cm"].to_numpy(float)}
        for year, code in ((2100, "21"), (2200, "22"), (2300, "23")):
            knots[year] = sampled[f"{v}_{code}"].to_numpy(float) * MM_TO_CM + shifts[v]
        xx = np.array([x2023, 2050.0, 2100.0, 2200.0, 2300.0])
        yy = np.vstack(
            [np.full(n, y2023), *[knots[y] for y in (2050, 2100, 2200, 2300)]]
        )
        spline = TensionSpline(xx, yy, p=tension_p, left_slope=rate2023)
        full[:, future] = spline.value(YEARS[future].astype(float)).T.astype(np.float32)
        for year, want in knots.items():
            actual = full[:, year - 2000].astype(float)
            diff = np.abs(actual - want)
            validations.append(
                {
                    "component": v,
                    "year": year,
                    "max_abs_difference_cm": float(diff.max()),
                    "allclose": bool(np.allclose(actual, want, rtol=1e-6, atol=1e-5)),
                }
            )
        for yr in (2050, 2100, 2200):
            diff = np.abs(spline.rate([yr + 1e-5])[0] - spline.rate([yr - 1e-5])[0])
            validations.append(
                {
                    "component": v,
                    "year": f"rate_{yr}",
                    "max_abs_difference_cm": float(diff.max()),
                    "allclose": bool(diff.max() < 1e-3),
                }
            )
        output[v] = full
    output["AIS"] = (output["E"] + output["W"]).astype(np.float32)
    output["SLR"] = (output["AIS"] + output["G"]).astype(np.float32)
    validation = pd.DataFrame(validations)
    if not validation["allclose"].all():
        raise RuntimeError(
            f"Failed knot/rate validations: {validation.loc[~validation.allclose].to_string(index=False)}"
        )
    return output, validation
