"""
03_laboratory_changes.py

Reproduces Table 1 (laboratory markers, baseline vs 2 weeks) and the
corresponding sentences in the Results and Discussion.

Writes output/03_table1.csv and output/03_laboratory.txt
"""
import os
import numpy as np
import pandas as pd
from common import load, paired_wilcoxon, fmt_p, round_half_up, OUT

df = load()
lines = []


def out(s=''):
    print(s)
    lines.append(s)


MARKERS = [
    ('WBC', 'wbc', 'x10^9/L'),
    ('HGB', 'haemoglobin', 'g/L'),
    ('PLT', 'platelet', 'x10^9/L'),
    ('ALT', 'alt', 'U/L'),
    ('AST', 'ast', 'U/L'),
    ('LDH', 'ldh', 'U/L'),
    ('TB', 'total_bilirubin', 'umol/L'),
    ('TG', 'triglyceride', 'mmol/L'),
    ('Fib', 'fibrinogen', 'g/L'),
    ('Ferritin', 'ferritin', 'ug/L'),
    ('sCD25', 'scd25', 'pg/mL'),
]

rows = []
out('=' * 78)
out('TABLE 1. Laboratory markers, baseline to 2 weeks (Wilcoxon signed-rank, asymptotic)')
out('=' * 78)
out(f'{"Marker":<10}{"Baseline":>12}{"2 weeks":>12}{"n":>5}{"P":>12}')
for label, stem, unit in MARKERS:
    x = df[f'{stem}_baseline']
    y = df[f'{stem}_week2']
    m = x.notna() & y.notna()
    n, p = paired_wilcoxon(x[m], y[m], exact=False)
    rows.append({'marker': label, 'unit': unit,
                 'baseline_median': round(float(x[m].median()), 2),
                 'week2_median': round(float(y[m].median()), 2),
                 'n_pairs': n, 'p_value': (round(p, 4) if not np.isnan(p) else None)})
    out(f'{label:<10}{round_half_up(x[m].median()):>12}{round_half_up(y[m].median()):>12}{n:>5}{fmt_p(p):>12}')

pd.DataFrame(rows).to_csv(os.path.join(OUT, '03_table1.csv'), index=False, encoding='utf-8-sig')

# relative change — two different definitions, reported separately because the
# manuscript describes the reduction as a within-patient median
out('\n' + '-' * 78)
out('RELATIVE CHANGE: two definitions')
out('-' * 78)
for label, stem in (('Ferritin', 'ferritin'), ('sCD25', 'scd25')):
    x = df[f'{stem}_baseline']; y = df[f'{stem}_week2']
    m = x.notna() & y.notna()
    rel = ((x[m] - y[m]) / x[m] * 100)
    ratio = (x[m].median() - y[m].median()) / x[m].median() * 100
    out(f'{label}:')
    out(f'  median of within-patient percentage reductions = {rel.median():.1f}%')
    out(f'  reduction computed from the two medians         = {ratio:.1f}%')

with open(os.path.join(OUT, '03_laboratory.txt'), 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines))
print('\n-> output/03_table1.csv, output/03_laboratory.txt')
