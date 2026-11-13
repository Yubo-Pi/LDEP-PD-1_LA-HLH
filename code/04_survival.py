"""
04_survival.py

Reproduces the overall survival analysis (manuscript Results, "Survival";
Figure 4): Kaplan-Meier estimates with log-log confidence intervals, numbers at
risk, the labelled post hoc comparison by 2-week response, and the 4-month
landmark analyses.

Writes output/04_survival.txt
"""
import os
import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from common import load, km_loglog, OUT

df = load()
n = len(df)
t = df.survival_months.values
e = (df.survival_status == 'Dead').astype(int).values
lines = []


def out(s=''):
    print(s)
    lines.append(s)


out('=' * 74)
out('FOLLOW-UP')
out('=' * 74)
out(f'Median follow-up (median of individual follow-up times): {np.median(t):.1f} months'
    f' (range {t.min():.0f}-{t.max():.0f})')
kmr = KaplanMeierFitter().fit(t, 1 - e)
out(f'Median follow-up by reverse Kaplan-Meier: {kmr.median_survival_time_:.1f} months')
out(f'Deaths: {int(e.sum())}/{n} ({e.sum() / n * 100:.1f}%)')
out(f'Deaths within 4 months: {int(((t <= 4) & (e == 1)).sum())}')

out('\n' + '=' * 74)
out('TABLE 4A / FIGURE 4A — OVERALL SURVIVAL, ALL PATIENTS')
out('=' * 74)
out(f'{"Time":>8}{"OS":>10}{"95% CI":>18}{"at risk":>10}')
for tm in [1, 2, 3, 4, 6, 12, 24, 36]:
    S, lo, hi, risk = km_loglog(t, e, tm)
    out(f'{tm:>6} mo{S * 100:>9.1f}%{f"{lo*100:.1f}-{hi*100:.1f}":>18}{risk:>10}')

out('\n' + '=' * 74)
out('FIGURE 4B — BY 2-WEEK RESPONSE (post hoc, unadjusted)')
out('=' * 74)
for lab, sel in (('HLH-CR', df.hlh_response_2w == 'CR'),
                 ('HLH-PR', df.hlh_response_2w == 'PR'),
                 ('NR', df.hlh_response_2w == 'NR'),
                 ('PR or NR', df.hlh_response_2w != 'CR')):
    tt, ee = t[sel.values], e[sel.values]
    S, lo, hi, _ = km_loglog(tt, ee, 12)
    out(f'  {lab:<10} n={sel.sum():>2}  deaths={int(ee.sum())}  12-month OS {S*100:>5.1f}%'
        f'  (95% CI {lo*100:.1f}-{hi*100:.1f})')

cr = (df.hlh_response_2w == 'CR').values
res = logrank_test(t[cr], t[~cr], e[cr], e[~cr])
out(f'\n  log-rank, CR vs PR or NR: P = {res.p_value:.4f} (unadjusted, exploratory)')

out('\n' + '=' * 74)
out('4-MONTH LANDMARK — RESPONSE')
out('=' * 74)
alive4 = t > 4
out(f'Patients alive beyond 4 months (landmark set): {alive4.sum()}')
out(f'  not in the landmark set: {n - alive4.sum()} '
    f'(died within 4 months: {int(((t <= 4) & (e == 1)).sum())}; '
    f'last contact before 4 months: {int(((t <= 4) & (e == 0)).sum())})')
after = alive4 & (e == 1)
out(f'Deaths after the landmark: {int(after.sum())} (pid {df.loc[after, "pid"].tolist()})')
for lab, sel in (('CR', cr & alive4), ('PR', (df.hlh_response_2w == 'PR').values & alive4),
                 ('NR', (df.hlh_response_2w == 'NR').values & alive4)):
    out(f'  {lab:<4} in landmark set: {int(sel.sum()):>2}; deaths after landmark: {int((sel & after).sum())}')

out('\n' + '=' * 74)
out('4-MONTH LANDMARK — TRANSPLANTATION')
out('=' * 74)
hs = (df.allogeneic_hsct == 'Yes').values
out(f'Allogeneic HSCT: {int(hs.sum())}; deaths {int(e[hs].sum())}; alive {int((1 - e[hs]).sum())}')
out(f'No HSCT       : {int((~hs).sum())}; deaths {int(e[~hs].sum())}; alive {int((1 - e[~hs]).sum())}')
out(f'\nIn the landmark set:')
out(f'  HSCT        : {int((hs & alive4).sum())}; deaths after landmark {int((hs & alive4 & after).sum())}')
out(f'  no HSCT     : {int((~hs & alive4).sum())}; deaths after landmark {int((~hs & alive4 & after).sum())}')

out('\n' + '=' * 74)
out('CAUSES OF DEATH')
out('=' * 74)
for _, r in df[df.survival_status == 'Dead'].sort_values('survival_months').iterrows():
    out(f'  pid {int(r.pid):>2}  {r.survival_months:>5.1f} mo  {r.hlh_response_2w:<3}  {r.cause_of_death}')

with open(os.path.join(OUT, '04_survival.txt'), 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines))
print('\n-> output/04_survival.txt')
