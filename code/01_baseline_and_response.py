"""
01_baseline_and_response.py

Reproduces:
  - baseline description (manuscript Results, "Patient characteristics")
  - primary response analysis and Figure 2 numbers

Writes output/01_baseline.txt and output/01_response.txt
"""
import os
import numpy as np
import pandas as pd
from common import load, wilson, OUT

df = load()
n = len(df)
lines = []


def out(s=''):
    print(s)
    lines.append(s)


out('=' * 72)
out('BASELINE CHARACTERISTICS (n = %d)' % n)
out('=' * 72)

male = (df.sex == 'M').sum()
out(f'Sex: {male} male, {n - male} female')
out(f'Age, median (range): {df.age_years.median():.1f} ({df.age_years.min():.0f}-{df.age_years.max():.0f}) years')
out(f'ECOG >= 2: {(df.ecog >= 2).sum()}/{df.ecog.notna().sum()}')

out('\nDisease status at enrolment:')
for k, v in df.disease_status.value_counts().items():
    out(f'  {k:<16} {v}')

out('\nHistological subtype (WHO-HAEM5):')
for k, v in df.subtype_WHO_HAEM5.value_counts().items():
    out(f'  {k:<52} {v}')
out(f'  {"TOTAL":<52} {len(df)}')

out('\nLineage:')
for k, v in df.lineage.value_counts().items():
    lo, hi = wilson(v, n)
    out(f'  {k:<22} {v:>2} ({v / n * 100:.1f}%; 95% CI {lo:.1f}-{hi:.1f})')

out('\nAnn Arbor stage:')
for s in ['I', 'II', 'III', 'IV']:
    v = (df.ann_arbor_stage == s).sum()
    out(f'  Stage {s:<3} {v:>2} ({v / n * 100:.1f}%)')
adv = df.ann_arbor_stage.isin(['III', 'IV']).sum()
out(f'  Stage III-IV {adv:>2} ({adv / n * 100:.1f}%)')
out(f'  Stage not recorded: {df.ann_arbor_stage.isna().sum() + (df.ann_arbor_stage == "").sum()}')

out(f'\nB symptoms: {(df.b_symptoms == "Yes").sum()} ({(df.b_symptoms == "Yes").sum() / n * 100:.1f}%)')

ipi = df.ipi.dropna()
out(f'\nIPI available for {len(ipi)}/{n}')
out(f'  median {ipi.median():.0f} (range {ipi.min():.0f}-{ipi.max():.0f})')
for k in range(1, 6):
    v = (ipi == k).sum()
    if v:
        out(f'  IPI {k}: {v} ({v / len(ipi) * 100:.1f}%)')
out(f'  IPI >= 3: {(ipi >= 3).sum()} ({(ipi >= 3).sum() / len(ipi) * 100:.1f}%)')

out(f'\nPINK recorded: {df.pink.notna().sum()} (value(s): {sorted(df.pink.dropna().tolist())})')
out(f'IPS recorded:  {df.ips.notna().sum()} (value(s): {sorted(df.ips.dropna().tolist())})')

out('\nPlasma EBV-DNA at baseline:')
tested = df.ebv_dna_baseline_status.isin(['detectable', 'undetectable']).sum()
det = (df.ebv_dna_baseline_status == 'detectable').sum()
out(f'  tested {tested}; detectable {det} ({det / tested * 100:.1f}% of those tested)')
bd = df.loc[df.ebv_dna_baseline_status == 'detectable', 'ebv_dna_baseline']
out(f'  median among detectable {bd.median():,.0f} copies/mL (range {bd.min():,.0f}-{bd.max():,.0f})')
allv = df.ebv_dna_baseline.dropna()
out(f'  median among all tested (including zeros) {allv.median():,.0f} (range {allv.min():,.0f}-{allv.max():,.0f})')

out('\nPBMC EBV-DNA at baseline:')
pbmc = df.ebv_dna_pbmc_baseline.dropna()
pdet = (pbmc > 0).sum()
out(f'  detectable {pdet}/{len(pbmc)} ({pdet / len(pbmc) * 100:.1f}%)')
out(f'  median among detectable {pbmc[pbmc > 0].median():,.0f} copies/mL')
out(f'  median among all {pbmc.median():,.0f} (range {pbmc.min():,.0f}-{pbmc.max():,.0f})')

out('\nTreatment delivery:')
for k, v in df.cycles_delivered.value_counts().items():
    out(f'  {k} cycle(s): {v}')

# ---------------- responder ----------------
out('\n' + '=' * 72)
out('PRIMARY RESPONSE (2 weeks)')
out('=' * 72)
cr = (df.hlh_response_2w == 'CR').sum()
pr = (df.hlh_response_2w == 'PR').sum()
nr = (df.hlh_response_2w == 'NR').sum()
orr = cr + pr
for lab, k in (('CR', cr), ('PR', pr), ('NR', nr), ('HLH-ORR', orr)):
    lo, hi = wilson(k, n)
    out(f'  {lab:<8} {k:>2}/{n} = {k / n * 100:5.1f}%  (95% CI {lo:.1f}-{hi:.1f})')

out('\n4-week response (only patients assessed):')
w4 = df[df.hlh_response_4w.notna() & (df.hlh_response_4w != '')]
out(f'  assessed {len(w4)}')
for k, v in w4.hlh_response_4w.value_counts().items():
    out(f'  {k}: {v}')
out(f'  not assessed (bridged after 2 weeks): {n - len(w4)}')

with open(os.path.join(OUT, '01_baseline.txt'), 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines))
print('\n-> output/01_baseline.txt')
