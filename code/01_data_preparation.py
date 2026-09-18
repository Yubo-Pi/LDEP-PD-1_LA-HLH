#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
01_data_preparation.py
----------------------
Builds the de-identified, analysis-ready patient-level dataset for:

    Pi Y, et al. Adding asparaginase and PD-1 blockade to the DEP regimen
    achieves rapid disease control and Epstein-Barr virus clearance in
    EBV-positive lymphoma-associated hemophagocytic lymphohistiocytosis.
    ClinicalTrials.gov NCT05775705

Input  : raw workbook 'L-DEP-PD-1治LAHS原始数据CHECK版.xls' (sheet 1)
Output : ../data/LDEP_PD1_LA-HLH_deidentified.csv
         ../outputs/data_preparation_log.txt

The raw workbook contains no direct identifiers (no name, medical-record
number, national ID, or calendar date columns). Calendar dates were never
recorded: 'weeks_dx_to_enrolment' is stored as an interval in weeks, not a date.

De-identification applied here:
  * subject identifiers replaced by sequential study numbers (1-24);
  * free-text fields retained as clinical descriptors only;
  * no transformation of any analysis variable.

Requires: xlrd (for the legacy .xls workbook)
"""
from __future__ import annotations

import os
import re
import sys
import csv

try:
    import xlrd
except ImportError:
    sys.exit("xlrd is required:  pip install xlrd==2.0.1")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
# Path to the source workbook. Pass it as the first CLI argument, or place it
# at data/raw/L-DEP-PD-1治LAHS原始数据CHECK版.xls
RAW = (sys.argv[1] if len(sys.argv) > 1 else
       os.path.join(ROOT, 'data', 'raw', 'L-DEP-PD-1治LAHS原始数据CHECK版.xls'))
OUT_CSV = os.path.join(ROOT, 'data', 'LDEP_PD1_LA-HLH_deidentified.csv')
LOG = os.path.join(ROOT, 'outputs', 'data_preparation_log.txt')

# ---------------------------------------------------------------- column map
# Raw sheet layout (0-based column index), 116 columns:
#   0-13   clinical / treatment
#   14-18  outcome (survival, response)
#   19-22  signs at diagnosis          (+57 -> week 2 ; +77 -> week 4)
#   23-36  baseline laboratory         (+57 -> week 2 ; +77 -> week 4)
#   37-73  baseline cytokine panel (37 analytes, baseline only)
#   74-75  baseline EBV-DNA plasma / PBMC (-> +20 week 2 ; +20 week 4)
#   76-95  week 2 ; 96-115 week 4
CYTOKINES = [
    'mip_1a', 'sdf_1a', 'il_27', 'il_1b', 'il_2', 'il_4', 'il_5', 'ip_10',
    'il_6', 'il_7', 'il_8', 'il_10', 'eotaxin', 'il_12p70', 'il_13', 'il_17a',
    'il_31', 'il_1ra', 'rantes', 'gm_csf', 'tnf_a', 'mip_1b', 'ifn_a',
    'mcp_1', 'il_9', 'tnf_b', 'gro_a', 'il_1a', 'il_23', 'il_15', 'il_21',
    'il_22', 'ifn_g', 'st2', 'cd163', 'cxcl9', 'il_18',
]
LAB = [  # (baseline col, name)
    (23, 'wbc'), (24, 'anc'), (25, 'hgb'), (26, 'plt'), (27, 'alt'),
    (28, 'ast'), (29, 'ldh'), (30, 'tbil'), (31, 'creatinine'),
    (32, 'triglycerides'), (33, 'ferritin'), (34, 'fibrinogen'),
    (35, 'nk_activity'), (36, 'scd25'),
]
SIGNS = [(19, 'fever'), (20, 'tmax'), (21, 'hemophagocytosis'), (22, 'splenomegaly')]

# ---- header guard -----------------------------------------------------------
# Column positions are verified against the header text before anything is read,
# so that any future edit of the source workbook fails loudly instead of silently
# shifting every variable by one.
HEADER_CHECKS = {
    0: '方案',
    1: '性别',
    2: '年龄',
    4: '原发病',
    7: 'ECOG',
    8: '入组前合并症',
    10: '给药疗程',
    11: '治疗不良反应',
    13: '序贯治疗',
    14: '生存情况',
    17: '疗效评估',
    19: '发热',
    23: 'WBC',
    36: 'sCD25',
    37: 'MIP-1a',
    73: 'IL-18',
    74: 'EBV-DNA',
    94: 'EBV-DNA',
    114: 'EBV-DNA',
}

MISSING = {'', '-', 'ND', 'NA', 'nd', '/', '未查', '未测', 'N/A'}

# Any calendar date appearing in a free-text field is redacted: the study
# collected intervals, not dates, and a date could narrow down a subject.
DATE_PAT = re.compile(
    r'\d{4}\s*[/\-.年]\s*\d{1,2}\s*[/\-.月]\s*\d{1,2}\s*日?'
    r'|\d{4}\s*年\s*\d{1,2}\s*月?'
    r'|\d{1,2}\s*月\s*\d{1,2}\s*日'
    r'|20\d{2}\s*[-/]\s*\d{1,2}')
REDACTIONS = []


def redact(text, col, subject):
    """Remove calendar dates from free-text fields."""
    if not text:
        return text
    new, k = DATE_PAT.subn('[date redacted]', text)
    if k:
        REDACTIONS.append((subject, col, k))
    return new


def to_num(v):
    """Numeric value, or None for any missing / non-numeric token."""
    if isinstance(v, float):
        return v
    if isinstance(v, int):
        return float(v)
    s = str(v).strip().replace(',', '').replace('U/m', '').strip()
    if s in MISSING:
        return None
    s = s.lstrip('<>=').strip()
    try:
        return float(s)
    except ValueError:
        return None


def to_txt(v):
    return str(v).strip()


def norm_regimen(v):
    """Raw values are inconsistently coded ('DEP+PD1' vs 'DEP+PD-1')."""
    s = to_txt(v).replace(' ', '').upper()
    if s.startswith('L-DEP'):
        return 'L-DEP+PD-1'
    if s.startswith('DEP'):
        return 'DEP+PD-1'
    return s or None


def main():
    if not os.path.exists(RAW):
        sys.exit(
            "Raw workbook not found.\n"
            "Place the source file at:\n  %s\n"
            "(The source workbook is not included in this public release: "
            "'L-DEP-PD-1治LAHS原始数据CHECK版.xls'.)" % RAW)

    book = xlrd.open_workbook(RAW)
    sh = book.sheet_by_index(0)

    # ---- verify the expected layout before reading a single value ------------
    ncol = sh.ncols
    if ncol != 116:
        sys.exit('Expected the 116-column source layout after removal of the '
                 'subsequent-therapy free-text column; received %d columns.' % ncol)
    for c, kw in HEADER_CHECKS.items():
        if c >= ncol:
            sys.exit('header check failed: workbook has %d columns, expected at least %d'
                     % (ncol, c + 1))
        hh = to_txt(sh.cell_value(0, c)) + '|' + to_txt(sh.cell_value(1, c))
        if kw not in hh:
            sys.exit('header check failed at column %d: expected %r, found %r\n'
                     'The source workbook layout may have changed - update the column '
                     'map in this script before proceeding.' % (c, kw, hh))
    print('header check passed: %d columns, all anchors matched' % ncol)

    rows = list(range(2, sh.nrows))          # data start on row index 2
    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    os.makedirs(os.path.dirname(LOG), exist_ok=True)

    fields = (['subject_id', 'regimen', 'sex', 'age', 'weeks_dx_to_enrolment',
               'underlying_disease', 'subtype_detail', 'disease_status', 'ecog',
               'comorbidity', 'comorbidity_detail', 'cycles',
               'adverse_event', 'adverse_event_detail', 'hsct',
               'survival_status', 'survival_months',
               'cause_of_death', 'lost_to_followup', 'response_2w', 'response_4w']
              + [f'{n}_{s}' for s in ('baseline',) for _, n in SIGNS]
              + [f'{n}_{s}' for s in ('baseline',) for _, n in LAB]
              + [f'cyt_{c}' for c in CYTOKINES]
              + ['ebv_plasma_baseline', 'ebv_pbmc_baseline']
              + [f'{n}_{s}' for s in ('week2', 'week4') for _, n in SIGNS]
              + [f'{n}_{s}' for s in ('week2', 'week4') for _, n in LAB]
              + ['ebv_plasma_week2', 'ebv_pbmc_week2',
                 'ebv_plasma_week4', 'ebv_pbmc_week4'])

    records = []
    for i, r in enumerate(rows, start=1):
        rec = {
            'subject_id': i,
            'regimen': norm_regimen(sh.cell_value(r, 0)),
            'sex': to_num(sh.cell_value(r, 1)),                 # 0 = male, 1 = female
            'age': to_num(sh.cell_value(r, 2)),
            'weeks_dx_to_enrolment': to_num(sh.cell_value(r, 3)),
            'underlying_disease': to_num(sh.cell_value(r, 4)),   # 0 = lymphoma, 1 = LPD
            'subtype_detail': redact(to_txt(sh.cell_value(r, 5)), 'subtype_detail', i),
            'disease_status': to_num(sh.cell_value(r, 6)),       # 0 new / 1 refractory / 2 relapsed
            'ecog': to_num(sh.cell_value(r, 7)),
            'comorbidity': to_num(sh.cell_value(r, 8)),
            'comorbidity_detail': redact(to_txt(sh.cell_value(r, 9)), 'comorbidity_detail', i),
            'cycles': to_num(sh.cell_value(r, 10)),              # 0 = one, 1 = two
            'adverse_event': to_num(sh.cell_value(r, 11)),
            'adverse_event_detail': redact(to_txt(sh.cell_value(r, 12)), 'adverse_event_detail', i),
            'hsct': to_num(sh.cell_value(r, 13)),
            # The source workbook no longer carries the free-text 'subsequent
            # lymphoma-directed therapy' column - the only field that held a
            # calendar date. No reported analysis depends on it.
            'survival_status': to_num(sh.cell_value(r, 14)),     # 0 alive / 1 dead / 2 lost
            'survival_months': to_num(sh.cell_value(r, 15)),
            'cause_of_death': redact(to_txt(sh.cell_value(r, 16)), 'cause_of_death', i),
            # The outcome note of five surviving subjects reads '失访' (lost to
            # follow-up); flagged explicitly rather than inferred downstream.
            'lost_to_followup': 1 if (to_num(sh.cell_value(r, 14)) == 0
                                      and '失访' in to_txt(sh.cell_value(r, 16))) else 0,
            'response_2w': to_num(sh.cell_value(r, 17)),         # 0 CR / 1 PR / 2 NR
            'response_4w': to_num(sh.cell_value(r, 18)),
        }
        for c, name in SIGNS:
            rec[f'{name}_baseline'] = to_num(sh.cell_value(r, c))
        for c, name in LAB:
            rec[f'{name}_baseline'] = to_num(sh.cell_value(r, c))
        for j, cy in enumerate(CYTOKINES):
            rec[f'cyt_{cy}'] = to_num(sh.cell_value(r, 37 + j))
        rec['ebv_plasma_baseline'] = to_num(sh.cell_value(r, 74))
        rec['ebv_pbmc_baseline'] = to_num(sh.cell_value(r, 75))
        for off, tag in ((57, 'week2'), (77, 'week4')):
            for c, name in SIGNS:
                rec[f'{name}_{tag}'] = to_num(sh.cell_value(r, c + off))
            for c, name in LAB:
                rec[f'{name}_{tag}'] = to_num(sh.cell_value(r, c + off))
        rec['ebv_plasma_week2'] = to_num(sh.cell_value(r, 94))
        rec['ebv_pbmc_week2'] = to_num(sh.cell_value(r, 95))
        rec['ebv_plasma_week4'] = to_num(sh.cell_value(r, 114))
        rec['ebv_pbmc_week4'] = to_num(sh.cell_value(r, 115))
        records.append(rec)

    with open(OUT_CSV, 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(records)

    lines = ['Data preparation log',
             '=' * 60,
             'raw workbook        : %s' % os.path.basename(RAW),
             'sheets with data    : 1 (sheets 2 and 3 are empty)',
             'raw dimensions      : %d patients x %d columns' % (len(rows), sh.ncols),
             'variables exported  : %d' % len(fields),
             'output              : %s' % os.path.relpath(OUT_CSV, ROOT),
             '',
             'No direct identifiers were present in the source workbook;',
             'no calendar dates were collected (intervals only).',
             '',
             'The free-text column "subsequent lymphoma-directed therapy" was',
             'removed from the source workbook by the authors, because it was the',
             'only field anywhere in the file that contained a calendar date.',
             'No reported analysis depends on it, and it is therefore absent from',
             'the released dataset.',
             '',
             'Date tokens redacted inside other exported free-text fields: %d' % len(REDACTIONS)]
    for subj, col, k in REDACTIONS:
        lines.append('  subject %d, %s: %d date token(s) replaced with [date redacted]'
                     % (subj, col, k))
    lines.append('')
    lines.append('Note: %d surviving subjects are recorded as lost to follow-up '
                 '(lost_to_followup = 1).' % sum(r['lost_to_followup'] for r in records))
    with open(LOG, 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(lines) + '\n')

    print('wrote %s  (%d patients x %d variables)' % (OUT_CSV, len(records), len(fields)))
    print('log   %s' % LOG)


if __name__ == '__main__':
    main()
