"""
05_tables_3_4.py

Reproduces the legacy subtype-level outcome table (manuscript Table 3 in
dataset v1.1, dropped from the v1.2 manuscript) and the revised patient-level
Table 3 from the shared dataset.

Note (v1.2): three patients previously coded as 'EBV-positive T/NK-cell LPD,
grade 2-3' are now 'systemic CAEBVD, T/NK-cell type (category A2-A3)', and
three patients previously coded as 'T-cell lymphoma, NOS' / 'EBV-positive
T-cell lymphoma, NOS' are now 'PTCL, NOS', per WHO-HAEM5 entities used in the
revised manuscript.

Writes output/05_table3.csv and output/05_table4.csv
"""
import os
import numpy as np
import pandas as pd
from common import load, km_loglog, OUT

df = load()
lines = []


def out(s=''):
    print(s)
    lines.append(s)


ORDER = [
    'ENKTL',
    'PTCL, NOS',
    'systemic CAEBVD, T/NK-cell type (category A2-A3)',
    'EBV-positive diffuse large B-cell lymphoma, NOS',
    'ALCL, ALK-negative',
    'Classical Hodgkin lymphoma',
]
SHORT = {
    'ENKTL': 'ENKTL',
    'PTCL, NOS': 'PTCL, NOS',
    'systemic CAEBVD, T/NK-cell type (category A2-A3)': 'systemic CAEBVD, T/NK-cell type (A2-A3)',
    'EBV-positive diffuse large B-cell lymphoma, NOS': 'EBV+ DLBCL, NOS',
    'ALCL, ALK-negative': 'ALCL, ALK-negative',
    'Classical Hodgkin lymphoma': 'Classical Hodgkin lymphoma',
}

rows = []


def summarise(sub):
    n = len(sub)
    cr = (sub.hlh_response_2w == 'CR').sum()
    pr = (sub.hlh_response_2w == 'PR').sum()
    det = sub.ebv_dna_baseline_status == 'detectable'
    n_det = det.sum()
    cleared = (det & sub.ebv_dna_week2.notna() & (sub.ebv_dna_week2 == 0)).sum()
    hsct = (sub.allogeneic_hsct == 'Yes').sum()
    dead = (sub.survival_status == 'Dead').sum()
    return n, cr, pr, f'{cleared}/{n_det}', hsct, dead


out('=' * 92)
out('TABLE 3. Response, plasma EBV-DNA, transplantation and survival by histological subgroup')
out('=' * 92)
out(f'{"Histology":<34}{"n":>4}{"CR":>4}{"PR":>4}{"EBV cleared":>14}{"HSCT":>6}{"Dead":>6}{"1-y OS (95% CI)":>20}')

for sub_name in ORDER:
    sub = df[df.subtype_WHO_HAEM5 == sub_name]
    if len(sub) == 0:
        continue
    n, cr, pr, ebv, hsct, dead = summarise(sub)
    if dead == 0:
        os_str = '-'
    else:
        S, lo, hi, _ = km_loglog(sub.survival_months.values,
                                 (sub.survival_status == 'Dead').astype(int).values, 12)
        os_str = f'{S * 100:.0f} ({lo * 100:.0f}-{hi * 100:.0f})'
    rows.append({'Histology': SHORT[sub_name], 'n': n, 'CR_2w': int(cr), 'PR_2w': int(pr),
                 'EBV_cleared': ebv, 'HSCT': int(hsct), 'Dead': int(dead), 'OS_1y': os_str})
    out(f'{SHORT[sub_name]:<34}{n:>4}{cr:>4}{pr:>4}{ebv:>14}{hsct:>6}{dead:>6}{os_str:>20}')

# subtotal: T- or NK-cell lineage
tnk = df[df.lineage == 'T- or NK-cell']
n, cr, pr, ebv, hsct, dead = summarise(tnk)
S, lo, hi, _ = km_loglog(tnk.survival_months.values,
                         (tnk.survival_status == 'Dead').astype(int).values, 12)
sub_os = f'{S * 100:.0f} ({lo * 100:.0f}-{hi * 100:.0f})'
out(f'{"T- or NK-cell neoplasms, subtotal":<34}{n:>4}{cr:>4}{pr:>4}{ebv:>14}{hsct:>6}{dead:>6}{sub_os:>20}')
rows.append({'Histology': 'T- or NK-cell neoplasms, subtotal', 'n': n, 'CR_2w': int(cr),
             'PR_2w': int(pr), 'EBV_cleared': ebv, 'HSCT': int(hsct), 'Dead': int(dead),
             'OS_1y': sub_os})

n, cr, pr, ebv, hsct, dead = summarise(df)
S, lo, hi, _ = km_loglog(df.survival_months.values,
                         (df.survival_status == 'Dead').astype(int).values, 12)
all_os = f'{S * 100:.0f} ({lo * 100:.0f}-{hi * 100:.0f})'
out(f'{"All patients":<34}{n:>4}{cr:>4}{pr:>4}{ebv:>14}{hsct:>6}{dead:>6}{all_os:>20}')
rows.append({'Histology': 'All patients', 'n': n, 'CR_2w': int(cr), 'PR_2w': int(pr),
             'EBV_cleared': ebv, 'HSCT': int(hsct), 'Dead': int(dead), 'OS_1y': all_os})

pd.DataFrame(rows).to_csv(os.path.join(OUT, '05_table3.csv'), index=False, encoding='utf-8-sig')

out('\n' + '=' * 92)
out('TABLE 4. Patient-level data')
out('=' * 92)
t4 = pd.DataFrame({
    'No.': df.pid.astype(int),
    'Age/Sex': df.age_years.astype(int).astype(str) + '/' + df.sex,
    'Subtype': df.subtype_WHO_HAEM5,
    'Status': df.disease_status,
    'ECOG': df.ecog,
    'Response 2w': df.hlh_response_2w,
    'Response 4w': df.hlh_response_4w.fillna(''),
    'EBV base': df.ebv_dna_baseline,
    'EBV 2w': df.ebv_dna_week2,
    'Ferritin': df.ferritin_baseline,
    'sCD25': df.scd25_baseline,
    'HSCT': df.allogeneic_hsct,
    'OS (mo)': df.survival_months,
    'Outcome': df.survival_status,
})
t4.to_csv(os.path.join(OUT, '05_table4.csv'), index=False, encoding='utf-8-sig')
out(t4.to_string(index=False))

with open(os.path.join(OUT, '05_tables_3_4.txt'), 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines))
print('\n-> output/05_table3.csv, output/05_table4.csv')
