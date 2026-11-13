"""
02_ebv.py

Reproduces the plasma / PBMC EBV-DNA analysis (manuscript Results, "Changes in
EBV-DNA"; Figure 3) and the three denominators that the manuscript reports
separately.

Writes output/02_ebv.txt
"""
import os
import numpy as np
import pandas as pd
from common import load, wilson, paired_wilcoxon, round_half_up, OUT

df = load()
n = len(df)
lines = []


def out(s=''):
    print(s)
    lines.append(s)


b = df.ebv_dna_baseline
w = df.ebv_dna_week2

det = df.ebv_dna_baseline_status == 'detectable'
n_det = det.sum()
has2w = det & df.ebv_dna_week2_status.isin(['detectable', 'undetectable'])
cleared = det & (df.ebv_dna_week2 == 0) & df.ebv_dna_week2.notna()
n_has2w = has2w.sum()
n_cleared = cleared.sum()

out('=' * 72)
out('PLASMA EBV-DNA AT 2 WEEKS')
out('=' * 72)
out(f'Baseline tested            : {df.ebv_dna_baseline_status.isin(["detectable","undetectable"]).sum()}')
out(f'Baseline detectable        : {n_det}')
out(f'  of these, 2-week sample  : {n_has2w}')
out(f'  of these, no 2-week sample: {n_det - n_has2w}  (pid {df.loc[det & ~has2w, "pid"].tolist()})')
out(f'Baseline undetectable      : {(df.ebv_dna_baseline_status == "undetectable").sum()}'
    f'  (pid {df.loc[df.ebv_dna_baseline_status == "undetectable", "pid"].tolist()})')
out(f'Baseline not tested        : {(df.ebv_dna_baseline_status == "not_tested").sum()}'
    f'  (pid {df.loc[df.ebv_dna_baseline_status == "not_tested", "pid"].tolist()})')

lo, hi = wilson(n_cleared, n_det)
out(f'\n(1) Primary analysis, denominator = all with detectable baseline EBV-DNA,')
out(f'    patients without a 2-week sample counted as not cleared:')
out(f'    {n_cleared}/{n_det} = {round_half_up(n_cleared / n_det * 100)}%  (95% CI {lo:.1f}-{hi:.1f})')

lo2, hi2 = wilson(n_cleared, n_has2w)
out(f'\n(2) Complete-case analysis, denominator = detectable at baseline AND')
out(f'    a 2-week sample available:')
out(f'    {n_cleared}/{n_has2w} = {round_half_up(n_cleared / n_has2w * 100)}%  (95% CI {lo2:.1f}-{hi2:.1f})')

out(f'\n(3) Not at risk of clearance (already undetectable at baseline): '
    f'{(df.ebv_dna_baseline_status == "undetectable").sum()}')

# paired comparison
pair = df[df.ebv_dna_baseline.notna() & df.ebv_dna_week2.notna()]
out(f'\nPAIRED COMPARISON')
out(f'  patients with paired baseline and 2-week plasma samples: {len(pair)}')
out(f'  median baseline {pair.ebv_dna_baseline.median():,.0f} copies/mL')
out(f'  median at 2 weeks {pair.ebv_dna_week2.median():,.0f} copies/mL')
npair, pex = paired_wilcoxon(pair.ebv_dna_baseline, pair.ebv_dna_week2, exact=True)
out(f'  exact Wilcoxon signed-rank test: n = {npair} pairs, P = {pex:.4f}')
out(f'  (pairs with a zero difference are dropped by the exact test)')

# among paired with detectable baseline
paird = pair[pair.ebv_dna_baseline > 0]
out(f'\n  paired AND detectable at baseline: {len(paird)}')
out(f'  of these, undetectable at 2 weeks: {(paird.ebv_dna_week2 == 0).sum()}'
    f' ({round_half_up((paird.ebv_dna_week2 == 0).sum() / len(paird) * 100)}%)')

# waterfall categories
out('\nWATERFALL CATEGORIES (Figure 3)')
for k, v in df.ebv_plasma_2w_category.value_counts().items():
    pids = df.loc[df.ebv_plasma_2w_category == k, 'pid'].tolist()
    out(f'  {k:<28} n = {v:<2} pid {pids}')

# 4-week
out('\n4-WEEK PLASMA EBV-DNA')
w4 = df[df.ebv_dna_week4.notna()]
out(f'  tested: {len(w4)}; undetectable: {(w4.ebv_dna_week4 == 0).sum()}')

# PBMC
out('\n' + '=' * 72)
out('PBMC EBV-DNA')
out('=' * 72)
pb = df.ebv_dna_pbmc_baseline
pw = df.ebv_dna_pbmc_week2
out(f'  baseline detectable: {(pb > 0).sum()}/{pb.notna().sum()} ({round_half_up((pb > 0).sum() / pb.notna().sum() * 100)}%)')
out(f'  baseline median (all tested, as in the manuscript): {pb.median():,.0f} copies/mL'
    f' (range {pb.min():,.0f}-{pb.max():,.0f})')
out(f'  baseline median (detectable only): {pb[pb > 0].median():,.0f} copies/mL'
    f' (range {pb[pb > 0].min():,.0f}-{pb[pb > 0].max():,.0f})')
pairs = df[pb.notna() & pw.notna()]
out(f'  paired samples: {len(pairs)}')
out(f'  median at 2 weeks: {pairs.ebv_dna_pbmc_week2.median():,.0f} copies/mL')
npb, ppb = paired_wilcoxon(pairs.ebv_dna_pbmc_baseline, pairs.ebv_dna_pbmc_week2, exact=False)
out(f'  Wilcoxon signed-rank (asymptotic): n = {npb}, P = {ppb:.3f}')
pdet = pairs[pairs.ebv_dna_pbmc_baseline > 0]
out(f'  detectable at baseline among pairs: {len(pdet)};'
    f' undetectable at 2 weeks: {(pdet.ebv_dna_pbmc_week2 == 0).sum()}'
    f' ({(pdet.ebv_dna_pbmc_week2 == 0).sum() / len(pdet) * 100:.1f}%)')
pw4 = df[df.ebv_dna_pbmc_week4.notna()]
out(f'  4-week: tested {len(pw4)}; median {pw4.ebv_dna_pbmc_week4.median():,.0f};'
    f' undetectable {(pw4.ebv_dna_pbmc_week4 == 0).sum()}')

with open(os.path.join(OUT, '02_ebv.txt'), 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines))
print('\n-> output/02_ebv.txt')
