"""
common.py — shared helpers for the L-DEP-PD-1 LA-HLH analysis.

All scripts read data/ldep_pd1_lahlh_patient_level.csv and write to output/.
Run from the package root:  python code/run_all.py
"""
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, 'data', 'ldep_pd1_lahlh_patient_level.csv')
OUT = os.path.join(ROOT, 'output')
os.makedirs(OUT, exist_ok=True)

Z = 1.959963985


def load():
    df = pd.read_csv(DATA, encoding='utf-8-sig')
    df['pid'] = df['pid'].astype(int)
    return df


def wilson(k, n, z=Z):
    """Wilson score 95% CI, returned as percentages."""
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return (max(0.0, c - h) * 100, min(1.0, c + h) * 100)


def pct(k, n):
    return f'{k}/{n} ({k / n * 100:.1f}%'


def km_loglog(time, event, t_eval):
    """
    Kaplan-Meier estimate with Greenwood variance and a log-log transformed CI.
    Returns (S, lower, upper, number_at_risk) as proportions.
    """
    t = np.asarray(time, float)
    e = np.asarray(event, int)
    o = np.argsort(t)
    t, e = t[o], e[o]
    at_risk = len(t)
    S, cum = 1.0, 0.0
    last = (1.0, 0.0)
    for u in np.unique(t):
        d = int(((t == u) & (e == 1)).sum())
        c = int(((t == u) & (e == 0)).sum())
        if at_risk > 0 and d > 0:
            S *= (1 - d / at_risk)
            if at_risk > d:
                cum += d / (at_risk * (at_risk - d))
        if u <= t_eval:
            last = (S, cum)
        at_risk -= (d + c)
        if u > t_eval and d == 0:
            break
    S, cum = last
    if 0 < S < 1 and cum > 0:
        se = np.sqrt(cum) / abs(np.log(S))
        lo = S ** np.exp(Z * se)
        hi = S ** np.exp(-Z * se)
    else:
        lo = hi = S
    n_at = int((t >= t_eval).sum())
    return S, lo, hi, n_at


def paired_wilcoxon(x, y, exact=False):
    from scipy import stats
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = ~(np.isnan(x) | np.isnan(y))
    x, y = x[m], y[m]
    if len(x) < 2:
        return len(x), np.nan
    d = x - y
    nz = d != 0
    if nz.sum() < 1:
        return len(x), np.nan
    method = 'exact' if exact else 'approx'
    p = stats.wilcoxon(x, y, method=method).pvalue
    return len(x), float(p)


def fmt_p(p):
    """Report P to three decimals, with '<0.001' only below 0.0005."""
    if np.isnan(p):
        return 'NA'
    if p < 0.0005:
        return '<0.001'
    return f'{p:.3f}'


def round_half_up(v, nd=1):
    """Decimal half-up rounding, so that 2.15 -> 2.2 and 1088.45 -> 1088.5."""
    from decimal import Decimal, ROUND_HALF_UP
    q = Decimal(1).scaleb(-nd)
    return str(Decimal(str(round(float(v), 6))).quantize(q, rounding=ROUND_HALF_UP))
