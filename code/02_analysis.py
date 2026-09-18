#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
02_analysis.py
--------------
Reproduces every quantitative result reported in:

    Pi Y, et al. Adding asparaginase and PD-1 blockade to the DEP regimen
    achieves rapid disease control and Epstein-Barr virus clearance in
    EBV-positive lymphoma-associated hemophagocytic lymphohistiocytosis.
    ClinicalTrials.gov NCT05775705

Input  : ../data/LDEP_PD1_LA-HLH_deidentified.csv
Output : ../outputs/*.csv  (one file per reported table / result block)
         ../outputs/reproducibility_check.txt

Statistical conventions, matching the manuscript:
  * continuous variables are medians (and ranges);
  * paired baseline-to-week-2 comparisons use the two-sided asymptotic
    Wilcoxon signed-rank test (as stated in the manuscript's Statistical
    analysis section);
  * the paired plasma EBV-DNA comparison excludes zero-difference pairs and
    uses the exact Wilcoxon signed-rank test (n = 16);
  * binary proportions are reported with Clopper-Pearson (exact) 95% CIs;
  * overall survival is estimated by Kaplan-Meier, median follow-up as the
    median of individual follow-up times;
  * the assay lower limit of detection for EBV-DNA is 500 copies/mL; values
    at or below it are non-detectable.

Requires: pandas, numpy, scipy
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon, binomtest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, 'data', 'LDEP_PD1_LA-HLH_deidentified.csv')
OUT = os.path.join(ROOT, 'outputs')
LOD = 500.0                      # EBV-DNA assay lower limit of detection
RESP = {0: 'CR', 1: 'PR', 2: 'NR'}

REPORT = []


def say(line=''):
    print(line)
    REPORT.append(line)


def cp_ci(k, n):
    """Clopper-Pearson exact 95% CI for a proportion."""
    ci = binomtest(k, n).proportion_ci(0.95, method='exact')
    return 100 * ci.low, 100 * ci.high


def km(times, events):
    """Kaplan-Meier survival function, deaths ordered before censoring at ties."""
    data = sorted(zip(times, events), key=lambda z: z[0])
    at_risk = len(data)
    s = 1.0
    out_t, out_s, out_n = [], [], []
    i = 0
    while i < len(data):
        t = data[i][0]
        d = c = 0
        j = i
        while j < len(data) and data[j][0] == t:
            if data[j][1] == 1:
                d += 1
            else:
                c += 1
            j += 1
        if d and at_risk:
            s *= (1 - d / at_risk)
        at_risk -= (d + c)
        out_t.append(t); out_s.append(s); out_n.append(at_risk)
        i = j
    return np.array(out_t), np.array(out_s), np.array(out_n)


def km_at(t, s, tp):
    idx = np.where(t <= tp)[0]
    return s[idx[-1]] if len(idx) else 1.0


def wilcox_pair(df, a, b, method='asymptotic'):
    """Paired two-sided Wilcoxon signed-rank test on complete pairs."""
    d = df[[a, b]].dropna()
    if len(d) < 3:
        return len(d), np.nan, np.nan, np.nan
    x, y = d[a].values, d[b].values
    p = wilcoxon(x, y, alternative='two-sided', method=method).pvalue
    return len(d), float(np.median(x)), float(np.median(y)), float(p)


def main():
    if not os.path.exists(DATA):
        sys.exit('Dataset not found: %s\nRun 01_data_preparation.py first.' % DATA)
    os.makedirs(OUT, exist_ok=True)
    df = pd.read_csv(DATA)
    n = len(df)

    say('=' * 74)
    say('Reproducibility check - L-DEP-PD-1 in EBV-positive LA-HLH (NCT05775705)')
    say('dataset: %s' % os.path.basename(DATA))
    say('n = %d subjects, %d variables' % (n, df.shape[1]))
    say('=' * 74)

    # ------------------------------------------------------------------ 1
    say('')
    say('1. Patient characteristics')
    say('-' * 74)
    male = int((df.sex == 0).sum()); female = int((df.sex == 1).sum())
    say('  sex: %d male / %d female' % (male, female))
    say('  age: median %.1f (range %.0f-%.0f)' % (
        df.age.median(), df.age.min(), df.age.max()))
    say('  underlying disease: lymphoma %d, LPD %d' % (
        int((df.underlying_disease == 0).sum()), int((df.underlying_disease == 1).sum())))
    ds = df.disease_status.value_counts().to_dict()
    say('  disease status: newly diagnosed %d, refractory %d, relapsed %d' % (
        ds.get(0.0, 0), ds.get(1.0, 0), ds.get(2.0, 0)))
    say('  ECOG >= 2: %d' % int((df.ecog >= 2).sum()))
    say('  full L-DEP+PD-1 regimen: %d ; asparaginase omitted (DEP+PD-1): %d' % (
        int((df.regimen == 'L-DEP+PD-1').sum()), int((df.regimen == 'DEP+PD-1').sum())))
    say('  cycles given: one %d, two %d' % (
        int((df.cycles == 0).sum()), int((df.cycles == 1).sum())))

    bp = df.ebv_plasma_baseline.dropna()
    tested = len(bp); detectable = int((bp > LOD).sum())
    say('  plasma EBV-DNA at baseline: tested %d, detectable %d (%.1f%%), '
        'median %.0f (range %.0f-%.0f)' % (
            tested, detectable, 100 * detectable / tested,
            bp.median(), bp.min(), bp.max()))

    # ------------------------------------------------------------------ 2
    say('')
    say('2. Response at week 2')
    say('-' * 74)
    r = df.response_2w.dropna()
    cr = int((r == 0).sum()); pr = int((r == 1).sum()); nr = int((r == 2).sum())
    orr = cr + pr
    say('  CR %d (%.1f%%, 95%% CI %.1f-%.1f)' % (cr, 100 * cr / n, *cp_ci(cr, n)))
    say('  PR %d (%.1f%%)' % (pr, 100 * pr / n))
    say('  NR %d (%.1f%%)' % (nr, 100 * nr / n))
    say('  ORR %d/%d (%.1f%%, 95%% CI %.1f-%.1f)  [manuscript 95.8%%, 78.9-99.9]'
        % (orr, n, 100 * orr / n, *cp_ci(orr, n)))
    r4 = df.response_4w.dropna()
    if len(r4):
        o4 = int((r4 <= 1).sum())
        say('  week 4 evaluable %d: CR %d, PR %d, ORR %.1f%%' % (
            len(r4), int((r4 == 0).sum()), int((r4 == 1).sum()), 100 * o4 / len(r4)))

    pd.DataFrame([dict(metric='ORR_week2', k=orr, n=n, pct=100 * orr / n, ci_low=cp_ci(orr, n)[0], ci_high=cp_ci(orr, n)[1]),
                  dict(metric='CR_week2', k=cr, n=n, pct=100 * cr / n, ci_low=cp_ci(cr, n)[0], ci_high=cp_ci(cr, n)[1]),
                  dict(metric='PR_week2', k=pr, n=n, pct=100 * pr / n, ci_low=np.nan, ci_high=np.nan),
                  dict(metric='NR_week2', k=nr, n=n, pct=100 * nr / n, ci_low=np.nan, ci_high=np.nan)]
                 ).to_csv(os.path.join(OUT, 'result_response_week2.csv'), index=False)

    # ------------------------------------------------------------------ 3
    say('')
    say('3. EBV-DNA kinetics')
    say('-' * 74)
    for tag, label in (('plasma', 'plasma'), ('pbmc', 'PBMC')):
        a, b = f'ebv_{tag}_baseline', f'ebv_{tag}_week2'
        d = df[[a, b]].dropna()
        say('  [%s] paired n = %d ; median %.0f -> %.0f' % (
            label, len(d), d[a].median(), d[b].median()))
    d = df[['ebv_plasma_baseline', 'ebv_plasma_week2']].dropna()
    nz = d[d.ebv_plasma_baseline != d.ebv_plasma_week2]
    p_exact = wilcoxon(nz.ebv_plasma_baseline, nz.ebv_plasma_week2,
                       alternative='two-sided', method='exact').pvalue
    say('  [plasma] zero-difference pairs excluded -> n = %d ; exact Wilcoxon P = %.4f'
        '  [manuscript 0.009]' % (len(nz), p_exact))
    base_pos = d[d.ebv_plasma_baseline > LOD]
    cleared = base_pos[base_pos.ebv_plasma_week2 <= LOD]
    say('  [plasma] detectable at baseline and paired n = %d ; cleared at week 2 '
        '%d (%.1f%%, 95%% CI %.1f-%.1f)  [manuscript 9/16, 56.3%%]' % (
            len(base_pos), len(cleared), 100 * len(cleared) / len(base_pos),
            *cp_ci(len(cleared), len(base_pos))))
    neg = d[d.ebv_plasma_week2 <= LOD]
    say('  [plasma] non-detectable at week 2 %d/%d (%.1f%%)  [manuscript 11/18, 61.1%%]'
        % (len(neg), len(d), 100 * len(neg) / len(d)))
    d4 = df[['ebv_plasma_baseline', 'ebv_plasma_week4']].dropna()
    d4p = d4[d4.ebv_plasma_baseline > LOD]
    if len(d4p):
        say('  [plasma] week 4: evaluable %d ; cleared %d (%.1f%%)  [manuscript 5/6]' % (
            len(d4p), int((d4p.ebv_plasma_week4 <= LOD).sum()),
            100 * (d4p.ebv_plasma_week4 <= LOD).sum() / len(d4p)))

    dp = df[['ebv_pbmc_baseline', 'ebv_pbmc_week2']].dropna()
    p_pb = wilcoxon(dp.ebv_pbmc_baseline, dp.ebv_pbmc_week2,
                    alternative='two-sided', method='asymptotic').pvalue
    say('  [PBMC] paired n = %d ; median %.0f -> %.0f ; asymptotic Wilcoxon P = %.3f'
        '  [manuscript 0.334]' % (len(dp), dp.ebv_pbmc_baseline.median(),
                                  dp.ebv_pbmc_week2.median(), p_pb))
    pbp = dp[dp.ebv_pbmc_baseline > LOD]
    say('  [PBMC] detectable at baseline and paired n = %d ; cleared %d (%.1f%%)'
        '  [manuscript 9/18, 50%%]' % (
            len(pbp), int((pbp.ebv_pbmc_week2 <= LOD).sum()),
            100 * (pbp.ebv_pbmc_week2 <= LOD).sum() / len(pbp)))
    d4 = df[['ebv_pbmc_week4']].dropna()
    if len(d4):
        say('  [PBMC] week 4: evaluable %d ; median %.0f ; non-detectable %d'
            '  [manuscript 6, 2350, 3]' % (
                len(d4), d4.ebv_pbmc_week4.median(),
                int((d4.ebv_pbmc_week4 <= LOD).sum())))

    rows = []
    for tag, label, method in (('plasma', 'plasma', 'exact'), ('pbmc', 'PBMC', 'asymptotic')):
        a, b = f'ebv_{tag}_baseline', f'ebv_{tag}_week2'
        dd = df[[a, b]].dropna()
        rows.append(dict(compartment=label, n_paired=len(dd),
                         median_baseline=dd[a].median(), median_week2=dd[b].median(),
                         p_value=None))
    rows[0]['p_value'] = p_exact
    rows[1]['p_value'] = p_pb
    pd.DataFrame(rows).to_csv(os.path.join(OUT, 'result_ebv_kinetics.csv'), index=False)

    # ------------------------------------------------------------------ 4
    say('')
    say('4. Laboratory markers, baseline to week 2 (Table 1 of the manuscript)')
    say('-' * 74)
    say('  %-14s %12s %12s %5s %10s  %s' % (
        'marker', 'baseline', 'week2', 'n', 'P', 'manuscript'))
    MANU = {
        'wbc': ('1.6', '4.5', 24, '0.001'), 'hgb': ('90.0', '94.0', 24, '0.189'),
        'plt': ('47.5', '85.5', 24, '0.034'), 'alt': ('58.0', '35.0', 23, '0.001'),
        'ast': ('86.0', '26.7', 23, '0.001'), 'ldh': ('599.0', '322.0', 21, '0.004'),
        'tbil': ('20.9', '18.8', 24, '0.097'), 'triglycerides': ('2.2', '1.6', 15, '0.427'),
        'fibrinogen': ('1.5', '1.7', 24, '0.819'), 'ferritin': ('3894.2', '1088.5', 24, '<0.001'),
        'scd25': ('17391.0', '3406.0', 15, '0.001'),
    }
    lab_rows = []
    for m, (mb, m2, mn, mp) in MANU.items():
        npair, b, w, p = wilcox_pair(df, f'{m}_baseline', f'{m}_week2')
        say('  %-14s %12.1f %12.1f %5d %10.4f  %s->%s n=%d P=%s' % (
            m, b, w, npair, p, mb, m2, mn, mp))
        lab_rows.append(dict(marker=m, baseline_median=round(b, 1), week2_median=round(w, 1),
                             n_paired=npair, p_value=round(p, 4),
                             manuscript_baseline=mb, manuscript_week2=m2,
                             manuscript_n=mn, manuscript_p=mp))
    pd.DataFrame(lab_rows).to_csv(os.path.join(OUT, 'table1_laboratory.csv'), index=False)

    # ------------------------------------------------------------------ 5
    say('')
    say('5. Adverse events (Table 2 of the manuscript)')
    say('-' * 74)
    ae = int((df.adverse_event == 1).sum())
    say('  any adverse event: %d/%d (%.1f%%)  [manuscript 19/24, 79.2%%]'
        % (ae, n, 100 * ae / n))
    grades = {}
    for g in ('IV', 'III', 'II'):
        grades[g] = int(df.adverse_event_detail.fillna('').str.contains(
            g, regex=False).sum())
    # grade per subject = highest grade mentioned
    def top_grade(s):
        s = str(s)
        for g in ('IV', 'III', 'II'):
            if g in s:
                return g
        return 'none'
    tg = df.adverse_event_detail.fillna('').map(top_grade)
    say('  highest grade per subject: IV %d (%.1f%%), III %d (%.1f%%), II %d (%.1f%%)'
        '  [manuscript 12/5/2]' % (
            (tg == 'IV').sum(), 100 * (tg == 'IV').sum() / n,
            (tg == 'III').sum(), 100 * (tg == 'III').sum() / n,
            (tg == 'II').sum(), 100 * (tg == 'II').sum() / n))
    say('  no adverse event: %d' % int((tg == 'none').sum()))
    txt = df.adverse_event_detail.fillna('')
    for kw, label in (('肺炎', 'pneumonia'), ('巨细胞', 'CMV'), ('心房颤动', 'atrial fibrillation'),
                      ('室性早搏', 'ventricular premature beats')):
        say('  %-28s %d  [manuscript %s]' % (label, int(txt.str.contains(kw).sum()),
                                             {'pneumonia': 2, 'CMV': 2,
                                              'atrial fibrillation': 1,
                                              'ventricular premature beats': 1}[label]))
    pd.DataFrame([dict(category='any adverse event', n=ae),
                  dict(category='grade IV', n=int((tg == 'IV').sum())),
                  dict(category='grade III', n=int((tg == 'III').sum())),
                  dict(category='grade II', n=int((tg == 'II').sum())),
                  dict(category='pneumonia', n=int(txt.str.contains('肺炎').sum())),
                  dict(category='CMV', n=int(txt.str.contains('巨细胞').sum())),
                  dict(category='atrial fibrillation', n=int(txt.str.contains('心房颤动').sum())),
                  dict(category='ventricular premature beats', n=int(txt.str.contains('室性早搏').sum())),
                  dict(category='asparaginase omitted', n=int((df.regimen == 'DEP+PD-1').sum()))]
                 ).to_csv(os.path.join(OUT, 'table2_adverse_events.csv'), index=False)

    # ------------------------------------------------------------------ 6
    say('')
    say('6. Survival')
    say('-' * 74)
    t = df.survival_months.values.astype(float)
    e = (df.survival_status == 1).astype(int).values
    say('  deaths %d/%d (%.1f%%)  [manuscript 9/24, 37.5%%]' % (e.sum(), n, 100 * e.sum() / n))
    say('  follow-up: median %.1f months (range %.0f-%.0f)  [manuscript 8.5, 1-36]'
        % (np.median(t), t.min(), t.max()))
    kt, ks, _ = km(t, e)
    for tp, claim in ((3, '73.9%'), (6, '64.7%'), (12, '58.2%')):
        say('  OS at %2d months = %.1f%%  [manuscript %s]' % (tp, 100 * km_at(kt, ks, tp), claim))
    say('  deaths within 4 months: %d  [manuscript 8]' % int(((t <= 4) & (e == 1)).sum()))
    say('  surviving subjects lost to follow-up: %d' % int(df.lost_to_followup.sum()))
    say('')
    say('  cause of death (free text):')
    for _, row in df[df.survival_status == 1].iterrows():
        say('    subject %2d (OS %4.1f mo): %s' % (
            row.subject_id, row.survival_months, str(row.cause_of_death)[:70]))
    pd.DataFrame(dict(time=kt, survival=ks)).to_csv(
        os.path.join(OUT, 'result_survival_km.csv'), index=False)

    # ------------------------------------------------------------------ 7
    say('')
    say('7. Patient-level data (Table 3 of the manuscript)')
    say('-' * 74)
    t3 = df[['subject_id', 'age', 'sex', 'subtype_detail', 'disease_status', 'ecog',
             'response_2w', 'response_4w', 'ebv_plasma_baseline', 'ebv_plasma_week2',
             'ferritin_baseline', 'scd25_baseline', 'hsct', 'survival_months',
             'survival_status']].copy()
    t3['sex'] = t3.sex.map({0: 'M', 1: 'F'})
    t3['response_2w'] = t3.response_2w.map(RESP)
    t3['response_4w'] = t3.response_4w.map(RESP).fillna('-')
    t3.columns = ['No', 'Age', 'Sex', 'Subtype', 'Status', 'ECOG', 'Response_2w',
                  'Response_4w', 'EBV_DNA_baseline', 'EBV_DNA_week2', 'Ferritin',
                  'sCD25', 'HSCT', 'OS_months', 'Outcome']
    t3.to_csv(os.path.join(OUT, 'table3_patient_level.csv'), index=False)
    say('  written: table3_patient_level.csv (%d rows)' % len(t3))

    # ------------------------------------------------------------------ 8
    say('')
    say('8. Cross-checks against the manuscript narrative')
    say('-' * 74)
    d = df[['ebv_plasma_baseline', 'ebv_plasma_week2']].dropna()
    cr_ids = df.loc[df.response_2w == 0, 'subject_id'].tolist()
    ok, cleared_cr = 0, 0
    for sid in cr_ids:
        row = d.loc[df.subject_id == sid]
        if row.empty:
            continue
        b = float(row.ebv_plasma_baseline.iloc[0]); w = float(row.ebv_plasma_week2.iloc[0])
        if b > 0 and np.log10(max(b, LOD) / max(w, LOD)) >= 1:
            ok += 1
        if w <= LOD:
            cleared_cr += 1
    say('  CR at week 2: n = %d' % len(cr_ids))
    say('  ...of whom >= 1 log10 fall in plasma EBV-DNA: %d  [manuscript: all seven]' % ok)
    say('  ...of whom plasma EBV-DNA non-detectable at week 2: %d  '
        '[manuscript: four]' % cleared_cr)
    nr_id = df.loc[df.response_2w == 2, 'subject_id'].tolist()
    if nr_id:
        row = df[df.subject_id == nr_id[0]]
        say('  the single NR (subject %d): plasma EBV-DNA %.0f -> %.0f copies/mL'
            '  [manuscript: 270,000 -> 520,000]' % (
                nr_id[0], row.ebv_plasma_baseline.iloc[0], row.ebv_plasma_week2.iloc[0]))

    with open(os.path.join(OUT, 'reproducibility_check.txt'), 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(REPORT) + '\n')
    print('\nAll result files written to %s' % OUT)


if __name__ == '__main__':
    main()
