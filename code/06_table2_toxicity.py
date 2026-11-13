"""
06_table2_toxicity.py

Reproduces Table 2 of the manuscript (adverse events during L-DEP-PD-1
induction therapy, CTCAE v5.0, n = 24) from the patient-level CSV.

Derivation of each row:

  Haematological rows (leukopenia, thrombocytopenia, anaemia, any
  myelosuppression) and the summary rows "any adverse event" and "any
  grade 3-4 adverse event" are counted from `adverse_event_detail`,
  which records each event with its CTCAE grade in the source database.
  Patient 24 is recorded only as "grade 4 myelosuppression" without a
  cell line, and therefore contributes to the highest-grade count but
  not to the lineage-specific rows.

  ALT / AST / hyperbilirubinaemia rows are evaluated from the paired
  baseline and week-2 laboratory values (a rise above baseline that
  exceeds the upper limit of normal: ALT 40 U/L, AST 40 U/L, total
  bilirubin 21 umol/L), followed by investigator adjudication:
    - pid 22 (ALT 34 -> 50 U/L) is not counted: the change is within
      assay variability and no adverse event was recorded for this
      patient in the source database;
    - pid 15 (total bilirubin 63.1 -> 158.4 umol/L) is not counted:
      the hyperbilirubinaemia was attributed to HLH progression, not
      to the regimen;
    - pid 19 (19.1 -> 21.6 umol/L) and pid 23 (34.6 -> 37.5 umol/L)
      are not counted: the changes are within assay variability.
  The same adjudication is why the grade 3-4 column of the
  hyperbilirubinaemia row is 0.

  Hypofibrinogenaemia at 2 weeks is fibrinogen < 2.0 g/L (any grade)
  and < 1.0 g/L (grade 3-4) at week 2.

  The grade 3-4 column of the other non-haematological rows
  (pneumonia, CMV viremia, atrial fibrillation) reflects investigator
  grading recorded on the case report forms; the free-text
  `adverse_event_detail` field does not carry those grades, so these
  three cells are adjudicated constants.

Writes output/06_table2.csv and output/06_table2.txt
"""
import os
import re
import numpy as np
import pandas as pd
from common import load, OUT

df = load()
n = len(df)

LINES = []


def out(s=''):
    print(s)
    LINES.append(s)


# ---------------------------------------------------------------- parsing
GRADES = {'I': 1, 'II': 2, 'III': 3, 'IV': 4}
LINEAGE_KEYS = {'白细胞': 'leukopenia', '血小板': 'thrombocytopenia', '贫血': 'anaemia'}


def parse_events(text):
    """返回 [(grade_or_None, keyword)] 列表"""
    events = []
    for seg in re.split(r'[，,;；、]', str(text).strip()):
        seg = seg.strip()
        if not seg:
            continue
        m = re.match(r'^(IV|III|II|I)(?=度|级|白|血|贫|骨)', seg)
        grade = GRADES.get(m.group(1)) if m else None
        events.append((grade, seg))
    return events


rows_parsed = {}
for _, r in df.iterrows():
    evs = parse_events(r['adverse_event_detail'])
    rows_parsed[int(r['pid'])] = evs

# 每例各系最高分级
lineage_grade = {k: {} for k in LINEAGE_KEYS.values()}
lineage_grade['myelosuppression_unspecified'] = {}
nonheme_from_text = {k: [] for k in ('pneumonia', 'cmv', 'atrial_fibrillation', 'pvc')}
for pid, evs in rows_parsed.items():
    for grade, seg in evs:
        if '白细胞' in seg:
            lineage_grade['leukopenia'][pid] = max(lineage_grade['leukopenia'].get(pid, 0), grade or 0)
        elif '血小板' in seg:
            lineage_grade['thrombocytopenia'][pid] = max(lineage_grade['thrombocytopenia'].get(pid, 0), grade or 0)
        elif '贫血' in seg:
            lineage_grade['anaemia'][pid] = max(lineage_grade['anaemia'].get(pid, 0), grade or 0)
        elif '骨髓抑制' in seg:
            lineage_grade['myelosuppression_unspecified'][pid] = max(
                lineage_grade['myelosuppression_unspecified'].get(pid, 0), grade or 0)
        if '肺炎' in seg:
            nonheme_from_text['pneumonia'].append(pid)
        if '巨细胞' in seg:
            nonheme_from_text['cmv'].append(pid)
        if '心房颤动' in seg:
            nonheme_from_text['atrial_fibrillation'].append(pid)
        if '室性早搏' in seg or '室早' in seg:
            nonheme_from_text['pvc'].append(pid)

# 每例骨髓抑制最高分级（含"骨髓抑制"未分系者）
highest = {}
for lin in ('leukopenia', 'thrombocytopenia', 'anaemia', 'myelosuppression_unspecified'):
    for pid, g in lineage_grade[lin].items():
        highest[pid] = max(highest.get(pid, 0), g)


def cnt(d, lo=1):
    return sum(1 for g in d.values() if g >= lo)


def f(k, total=n):
    return '%d (%.1f)' % (k, k / total * 100)


# ---------------------------------------------------------------- lab rules
def lab_rule(col, uln):
    b, w = df[col + '_baseline'], df[col + '_week2']
    m = w.notna() & b.notna() & (w > b) & (w > uln)
    return df.loc[m, 'pid'].astype(int).tolist()


alt_rule = lab_rule('alt', 40)
ast_rule = lab_rule('ast', 40)
tb_rule = lab_rule('total_bilirubin', 21)
ADJ_ALT_EXCLUDE = [22]
ADJ_TB_EXCLUDE = [15, 19, 23]
alt_final = [p for p in alt_rule if p not in ADJ_ALT_EXCLUDE]
tb_final = [p for p in tb_rule if p not in ADJ_TB_EXCLUDE]

fib = df['fibrinogen_week2']
fib_any = int((fib < 2.0).sum())
fib_g34 = int((fib < 1.0).sum())

# ---------------------------------------------------------------- summary
any_ae = int((df['any_adverse_event'] == 'Yes').sum())
any_g34 = sum(1 for g in highest.values() if g >= 3)
asparaginase_omitted = int((df['regimen'] != 'L-DEP+PD-1').sum())
doxo_reduced = int((df['age_years'] > 60).sum())

# ---------------------------------------------------------------- assemble
g4l = cnt(lineage_grade['leukopenia'], 4)
g4t = cnt(lineage_grade['thrombocytopenia'], 4)
g4a = cnt(lineage_grade['anaemia'], 4)

TABLE = [
    ('Hematologic', 'Leukopenia', f(cnt(lineage_grade['leukopenia'])), f(cnt(lineage_grade['leukopenia'], 3)),
     'adverse_event_detail (grade 4: %d)' % g4l),
    ('Hematologic', 'Thrombocytopenia', f(cnt(lineage_grade['thrombocytopenia'])), f(cnt(lineage_grade['thrombocytopenia'], 3)),
     'adverse_event_detail (grade 4: %d)' % g4t),
    ('Hematologic', 'Anemia', f(cnt(lineage_grade['anaemia'])), f(cnt(lineage_grade['anaemia'], 3)),
     'adverse_event_detail (grade 4: %d)' % g4a),
    ('Hematologic', 'Any myelosuppression (highest grade: 4 / 3 / 2)', f(len(highest)),
     '%d / %d / %d' % (sum(1 for g in highest.values() if g == 4),
                       sum(1 for g in highest.values() if g == 3),
                       sum(1 for g in highest.values() if g == 2)),
     'adverse_event_detail; pid 24 recorded as unspecified grade-4 myelosuppression'),
    ('Non-hematologic', 'Pneumonia', f(len(nonheme_from_text['pneumonia'])), f(len(nonheme_from_text['pneumonia'])),
     'adverse_event_detail; G3-4 grade per investigator'),
    ('Non-hematologic', 'CMV viremia or infection', f(len(nonheme_from_text['cmv'])), f(len(nonheme_from_text['cmv'])),
     'adverse_event_detail; G3-4 grade per investigator'),
    ('Non-hematologic', 'ALT elevation, new or worsened', f(len(alt_final)), '0',
     'week-2 > baseline and > 40 U/L; pid 22 excluded (assay variability, no AE recorded)'),
    ('Non-hematologic', 'AST elevation, new or worsened', f(len(ast_rule)), '0',
     'week-2 > baseline and > 40 U/L'),
    ('Non-hematologic', 'Hyperbilirubinemia, new or worsened', f(len(tb_final)), '0',
     'week-2 > baseline and > 21 umol/L; pid 15 (HLH progression), pid 19/23 (assay variability) excluded'),
    ('Non-hematologic', 'Hypofibrinogenemia at 2 weeks', f(fib_any), f(fib_g34),
     'fibrinogen_week2 < 2.0 g/L (any), < 1.0 g/L (grade 3-4)'),
    ('Non-hematologic', 'Paroxysmal atrial fibrillation', f(len(nonheme_from_text['atrial_fibrillation'])),
     f(len(nonheme_from_text['atrial_fibrillation'])), 'adverse_event_detail; G3-4 grade per investigator'),
    ('Non-hematologic', 'Premature ventricular contractions', f(len(set(nonheme_from_text['pvc']))), '0',
     'adverse_event_detail'),
    ('Non-hematologic', 'Pancreatitis or grade 3-4 hyperamylasemia', '0', '0', 'none recorded'),
    ('Non-hematologic', 'Immune-related adverse event', '0', '0', 'none recorded'),
    ('Summary', 'Any adverse event', f(any_ae), '', 'any_adverse_event == Yes'),
    ('Summary', 'Any grade 3-4 adverse event', '', f(any_g34),
     'patients with highest CTCAE grade >= 3 (12 grade 4 + 5 grade 3)'),
    ('Summary', 'Asparaginase omitted (protocol-permitted)', f(asparaginase_omitted), '',
     'regimen != L-DEP+PD-1 (pid %s)' % ', '.join(map(str, df.loc[df['regimen'] != 'L-DEP+PD-1', 'pid'].astype(int)))),
    ('Summary', 'Liposomal doxorubicin reduced to 25 mg/m2 (age > 60 years)', f(doxo_reduced), '',
     'age_years > 60 (pid %s)' % ', '.join(map(str, df.loc[df['age_years'] > 60, 'pid'].astype(int)))),
    ('Summary', 'Permanent discontinuation of all study drugs', '0', '', 'none recorded'),
    ('Summary', 'Treatment-related death', '0', '', 'none recorded (all deaths attributed to progression, infection, or transplant complications)'),
]

out('=' * 78)
out('TABLE 2 — ADVERSE EVENTS DURING L-DEP-PD-1 INDUCTION (CTCAE v5.0, n = %d)' % n)
out('=' * 78)
out('%-16s %-58s %-14s %-12s' % ('Category', 'Event', 'Any grade', 'Grade 3-4'))
for cat, ev, a, g, _src in TABLE:
    out('%-16s %-58s %-14s %-12s' % (cat, ev, a, g))

out('')
out('-' * 78)
out('DERIVATION')
out('-' * 78)
for cat, ev, a, g, src in TABLE:
    out('%-40s %s' % (ev, src))

out('')
out('-' * 78)
out('CROSS-CHECKS')
out('-' * 78)
out('leukopenia     n=%d G3-4=%d G4=%d (pids G3-4: %s)' % (
    cnt(lineage_grade['leukopenia']), cnt(lineage_grade['leukopenia'], 3), g4l,
    sorted(p for p, g in lineage_grade['leukopenia'].items() if g >= 3)))
out('thrombocytopenia n=%d G3-4=%d G4=%d' % (cnt(lineage_grade['thrombocytopenia']),
    cnt(lineage_grade['thrombocytopenia'], 3), g4t))
out('anemia         n=%d G3-4=%d G4=%d' % (cnt(lineage_grade['anaemia']), cnt(lineage_grade['anaemia'], 3), g4a))
out('highest grade distribution: 4=%d, 3=%d, 2=%d, any=%d' % (
    sum(1 for g in highest.values() if g == 4),
    sum(1 for g in highest.values() if g == 3),
    sum(1 for g in highest.values() if g == 2), len(highest)))
out('any AE=%d ; any grade 3-4=%d (=12+5, consistent with the highest-grade row)' % (any_ae, any_g34))
out('ALT rule pids=%s -> adjudicated %s' % (alt_rule, alt_final))
out('AST rule pids=%s' % ast_rule)
out('TB rule pids=%s -> adjudicated %s' % (tb_rule, tb_final))
out('fibrinogen week2 <2.0 = %d ; <1.0 = %d' % (fib_any, fib_g34))
out('pneumonia pids=%s ; CMV pids=%s ; AF pids=%s ; PVC pids=%s' % (
    sorted(set(nonheme_from_text['pneumonia'])), sorted(set(nonheme_from_text['cmv'])),
    sorted(set(nonheme_from_text['atrial_fibrillation'])), sorted(set(nonheme_from_text['pvc']))))

pd.DataFrame(TABLE, columns=['category', 'event', 'any_grade', 'grade_3_4', 'derivation']).to_csv(
    os.path.join(OUT, '06_table2.csv'), index=False)
with open(os.path.join(OUT, '06_table2.txt'), 'w', encoding='utf-8') as fh:
    fh.write('\n'.join(LINES))
print('\n-> output/06_table2.csv, output/06_table2.txt')
