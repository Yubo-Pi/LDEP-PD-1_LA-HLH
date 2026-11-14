# De-identified patient-level dataset and analysis code for L-DEP-PD-1 induction therapy in EBV-positive lymphoma-associated HLH

**Version 1.2** · released 14 November 2026

This repository accompanies the manuscript:

> Pi Y, Wang J, Zhou H, Pan X, Wang Z. *Addition of asparaginase and PD-1 blockade to the DEP regimen was associated with rapid disease control and undetectable plasma EBV-DNA in EBV-positive lymphoma-associated hemophagocytic lymphohistiocytosis: a prospective, multicenter, single-arm phase 2 study.*

Trial registration: **NCT05775705**
Ethics approval: Beijing Friendship Hospital Ethics Committee, **No. 2023-P2-040-01**

---

## 1. What the study is

A prospective, multicentre, single-arm phase 2 study of L-DEP-PD-1 as induction therapy for adults with Epstein-Barr virus (EBV)-positive lymphoma-associated haemophagocytic lymphohistiocytosis (LA-HLH). Twenty-four patients were enrolled at three centres in China between August 2023 and July 2026. The regimen combined liposomal doxorubicin, etoposide, methylprednisolone, asparaginase and camrelizumab, given every 2 weeks for up to two cycles. The primary endpoint was the HLH overall response rate (HLH-ORR) at 2 weeks.

The study was not designed to demonstrate superiority over the DEP backbone, and the data here cannot support such a claim.

---

## 2. Contents

```
LDEP-PD-1_LA-HLH_v1.2/
├── README.md                        this file
├── LICENSE-data.txt                 CC BY 4.0 (applies to data/)
├── LICENSE-code.txt                 MIT (applies to code/)
├── data/
│   ├── ldep_pd1_lahlh_patient_level.csv    24 rows x 126 columns
│   └── data_dictionary.csv                 one row per variable
├── code/
│   ├── common.py                    shared helpers (Wilson CI, KM with log-log CI)
│   ├── 01_baseline_and_response.py  baseline table and primary response
│   ├── 02_ebv.py                    plasma / PBMC EBV-DNA analysis
│   ├── 03_laboratory_changes.py     Table 1 and relative changes
│   ├── 04_survival.py               Figure 4, landmark analyses
│   ├── 05_tables_3_4.py             legacy subtype table; patient-level table
│   ├── 06_table2_toxicity.py        Table 2 (adverse events)
│   ├── run_all.py                   runs everything in order
│   └── requirements.txt
└── output/                          created when the code is run
```

---

## 3. The dataset

### 3.1 Structure

One row per patient, 24 rows. `pid` runs from 1 to 24 and matches the patient numbers in Table 3 of the manuscript. Variables are named in `snake_case`; the suffix `_baseline`, `_week2` or `_week4` identifies the time point. All 126 variables are described in `data/data_dictionary.csv`.

### 3.2 De-identification

The dataset was de-identified before release:

- names, hospital numbers, dates of birth and exact admission dates were removed;
- no free-text field contains a date, an address or an identifier;
- `pid` is a study-specific sequential number and cannot be linked back to hospital records;
- the `comorbidity` and `cause_of_death` fields retain the original clinical wording in Chinese because the clinical detail is essential to interpretation; both were reviewed and contain no identifying information.

Dates survive only as an interval (`weeks_from_diagnosis_to_enrolment`) and as survival time in months.

### 3.3 Conventions and missing data

| Element | Convention |
|---|---|
| Empty cell | Not measured, not recorded, or not applicable — see the paired `*_status` column for EBV |
| `ebv_dna_*_status` | `detectable`, `undetectable`, `not_tested` (sample never taken) or `missing` (tested at another site, result unavailable) |
| Plasma EBV-DNA | copies/mL. `0` means below the limit of quantification of 500 copies/mL; it does **not** mean confirmed absence of virus |
| `survival_status` | `Alive` (10), `Dead` (9) or `Lost to follow-up` (5). The source database recorded only alive/dead; the five patients lost to follow-up are identified by the string `失访` in `cause_of_death` and are restored here, so that the dataset matches Table 3 and Figure 2 of the manuscript (10 alive / 5 lost / 9 dead) |
| `survival_months` | from the first dose of L-DEP-PD-1 to death or last contact; for patients lost to follow-up it is the time to last contact. Patients lost to follow-up are censored in the Kaplan–Meier analysis |
| `hlh_response_2w` | modified Marsh criteria, as defined in the manuscript Methods |
| `hlh_response_4w` | empty when the patient started lymphoma-directed therapy before the 4-week assessment |
| `gr` | granulocyte count, recorded in the source database but not used in the manuscript |

### 3.3.1 Coding of `regimen`

`regimen` has two values only:

- `L-DEP+PD-1` — the full regimen including asparaginase (20 patients)
- `DEP+PD-1` — L-DEP-PD-1 with asparaginase omitted, which the protocol permitted for organ intolerance (4 patients: pid 1, 11, 21, 22)

In version 1.1 the record for pid 1 was entered as `DEP+PD1`, without the hyphen. Any analysis that selected on `regimen != 'L-DEP+PD-1'` was unaffected, because the string differs from the full-regimen label either way, but a selection on the literal `'DEP+PD-1'` would have returned three patients instead of four. Version 1.2 corrects the entry; no other value in the dataset changed.

### 3.3.2 Coding of the grade 3–4 column of `06_table2.csv`

Most cells of `06_table2.csv` are counted mechanically, either from the free-text `adverse_event_detail` field or from paired laboratory values, and the `derivation` column names the rule used for each row. Four non-haematological rows are set by investigator adjudication rather than counted:

- **Pneumonia** — 2 patients, grade 3–4 in 0. The event is counted, but the source record carries no CTCAE grade for it.
- **CMV viraemia or infection** — 2 patients, grade 3–4 in 0. Same basis.
- **Paroxysmal atrial fibrillation** — 1 patient, grade 3–4 in 0. Same basis.
- **Hypofibrinogenaemia at 2 weeks** — 18 patients below 2.0 g/L, grade 3–4 in 0. One patient sits at 0.99 g/L, just under the 1.0 g/L CTCAE grade 3 threshold; the fall was not graded as an adverse event, so the row is reported as not reaching grade 3–4 rather than counted from the value.

These four rows are held in the `ADJUDICATED_GRADE_3_4` dictionary at the top of `code/06_table2_toxicity.py` and match Table 2 of the manuscript. Do not recompute them from `adverse_event_detail` alone: a plain substring search on that field returns more patients with pneumonia than the table reports, because the field also mentions pneumonia in the comorbidity and past-history text. The count in the table reflects investigator adjudication of which events were treatment-emergent.

Because every non-haematological row now reads 0, the "Any grade 3-4 adverse event" figure of 17 (70.8%) can look inconsistent with the table. It is not. That figure counts patients whose *highest* CTCAE grade was 3 or 4, from haematological grades alone (12 grade 4 + 5 grade 3), and it includes patient 24, recorded only as unspecified grade-4 myelosuppression. The non-haematological rows and this total are not meant to add up.

### 3.4 Histological classification

The source database recorded histology as free text, with several spellings for the same entity. `subtype_WHO_HAEM5` gives the entity assigned on re-review of every pathology report against the fifth edition of the WHO classification, and matches Table 3 of the manuscript. `lineage` groups these into T- or NK-cell versus B-cell or Hodgkin neoplasms:

- **T- or NK-cell (21 patients):** 13 ENKTL, 4 PTCL NOS, 3 systemic CAEBVD of T/NK-cell type (category A2–A3), 1 ALK-negative ALCL
- **B-cell or Hodgkin (3 patients):** 2 EBV-positive DLBCL NOS, 1 classical Hodgkin lymphoma

### 3.5 Ann Arbor stage and IPI

Both were available for all 24 patients. Stage was parsed from the free-text field: IV 17 (70.8%), III 3 (12.5%), II 1 (4.2%), I 3 (12.5%). The IPI had a median of 3 (range 1–5), and 19 patients (79.2%) scored 3 or higher. PINK was recorded for one ENKTL patient and the IPS for the single patient with classical Hodgkin lymphoma.

---

## 4. Reproducing the analysis

### 4.1 Requirements

Python 3.11 or later, with:

```
pandas>=2.0
numpy>=1.24
scipy>=1.11
lifelines>=0.27
```

Install with `pip install -r code/requirements.txt`. The analysis was run on Python 3.13.12.

### 4.2 Run

```bash
cd LDEP-PD-1_LA-HLH_v1.2
python code/run_all.py
```

Individual scripts can be run on their own; `common.py` must be importable, so run them from the package root or with `code/` on `PYTHONPATH`.

### 4.3 Output

Each script writes a plain-text log and, where relevant, a CSV to `output/`:

| Script | Output | Manuscript element |
|---|---|---|
| 01 | `01_baseline.txt` | Results, "Patient characteristics" and "Treatment response"; Figure 2 |
| 02 | `02_ebv.txt` | Results, "Changes in EBV-DNA"; Figure 3 |
| 03 | `03_table1.csv`, `03_laboratory.txt` | Table 1; Results, "Changes in laboratory markers" |
| 04 | `04_survival.txt` | Results, "Survival"; Figure 4A and 4B |
| 05 | `05_table3.csv`, `05_table4.csv` | Legacy subtype summary table (dropped from the revised manuscript) and patient-level table (manuscript Table 3) |
| 06 | `06_table2.csv`, `06_table2.txt` | Table 2 (adverse events) |

### 4.4 Statistical conventions

- **Proportions** — Wilson score 95% confidence intervals, no continuity correction.
- **Paired comparisons** — Wilcoxon signed-rank test. Asymptotic (two-sided) for the laboratory markers in Table 1. Exact for plasma EBV-DNA, because the number of non-zero differences is small; pairs with a zero difference are dropped by the exact test, which is why that comparison rests on 16 of the 18 paired samples.
- **Survival** — Kaplan–Meier with Greenwood variance and a log-log transformed 95% confidence interval. Numbers at risk are counted as patients with follow-up at or beyond each time point.
- **Group comparison** — log-rank test, unadjusted. The comparison by 2-week response is a labelled post hoc analysis and the study was not powered for it.
- **Landmark analysis** — restricted to patients with follow-up beyond 4 months, to address early mortality and immortal time bias. This set contains 14 patients: 8 died within 4 months and 2 were lost to follow-up before 4 months.
- **`05_table3.csv`, 1-year OS by lymphoma subtype (legacy output dropped from the revised manuscript)** — the estimate is shown for every subgroup in which at least one death occurred. It is left as `-` where no death occurred, because Kaplan–Meier returns exactly 100% in that situation and the number carries no information. For subgroups of one or two patients the estimate is printed as computed, but the confidence interval spans almost the whole range and the value should not be interpreted.
- **Rounding** — medians are reported to one decimal place; percentages to one decimal place except in the legacy subtype table (`05_table3.csv`), where they are given as whole numbers. Rounding is half-up (so 1088.45 → 1088.5 and 56.25% → 56.3%). The CSV outputs (`03_table1.csv`) carry unrounded medians and four-decimal P values; the companion text outputs apply the manuscript's rounding and reproduce the tables as printed.

---

## 5. Key results

These are the values produced by `code/run_all.py` from the CSV in `data/`.

**Baseline.** 17 male, 7 female; median age 36.5 years (range 18–70). 21 patients had a T- or NK-cell neoplasm and 3 a B-cell or Hodgkin neoplasm. Stage III–IV disease in 20 (83.3%); IPI ≥ 3 in 19 (79.2%).

**Response at 2 weeks.** HLH-ORR 23/24 = 95.8% (95% CI 79.8–99.3); CR 7/24 = 29.2% (14.9–49.2); PR 16/24 = 66.7%; NR 1/24 = 4.2%.

**Plasma EBV-DNA.** Detectable at baseline in 20 of the 23 patients tested. Undetectable at 2 weeks in 9 of 20 = 45.0% (25.8–65.8) when patients without a 2-week sample are counted as not cleared, and 9 of 16 = 56.3% (33.2–76.9) in a complete-case analysis. Median fell from 31,000 to 0 copies/mL among the 18 patients with paired samples (P = 0.009, exact Wilcoxon, 16 non-zero pairs).

**Laboratory markers.** Ferritin 3894.2 → 1088.5 µg/L (P < 0.001); sCD25 17,391 → 3,406 pg/mL (P = 0.001). The percentage reductions quoted in the manuscript (72.0% and 80.4%) are the reductions of the median values; the median of the within-patient reductions is 75.9% and 91.9% respectively, and `output/03_laboratory.txt` reports both.

**Survival.** Median follow-up 8.5 months (range 1–36). Nine of 24 patients died (37.5%), eight within 4 months; five (20.8%) were lost to follow-up after a median of 7 months (range 1.5–12), and the remaining ten were alive at the last contact. Overall survival 91.7% at 1 month, 83.1% at 2 months, 73.9% at 3 months, 64.7% at 4 months, 58.2% at 12 months (95% CI 34.2–76.1), unchanged thereafter. Twelve-month survival was 80.0% after a 2-week CR and 51.4% after a PR; the difference did not reach significance (log-rank P = 0.089).

**By histological subgroup.** The 21 T- or NK-cell neoplasms had an HLH-ORR of 20/21, CR 6, EBV clearance 8/17, 9 transplants, 9 deaths, and 1-year OS 53% (95% CI 28–72).

**Toxicity (Table 2).** All counts below are produced by `06_table2_toxicity.py`. Nineteen of 24 patients (79.2%) had at least one adverse event, and 17 (70.8%) had at least one grade 3-4 event (12 patients with a highest grade of 4 plus 5 with a highest grade of 3, consistent with the "highest grade 4 / 3 / 2 = 12 / 5 / 2" row). The haematological rows are counted from `adverse_event_detail`: leukopenia 18/24 (75.0%), grade 3-4 in 15 (62.5%), grade 4 in 9 (37.5%); thrombocytopenia 10/24 (41.7%), grade 3-4 in 10, grade 4 in 8 (33.3%); anaemia 5/24 (20.8%), all grade 3-4, grade 4 in 3. Patient 24 is recorded only as "grade 4 myelosuppression" without a cell line and contributes to the highest-grade count but not to the lineage-specific rows. Non-haematological events recorded in the same field: pneumonia 2, CMV viremia or infection 2, paroxysmal atrial fibrillation 1, premature ventricular contractions 1. ALT and AST elevations are evaluated from the paired laboratory values (week-2 above baseline and above the upper limit of normal), giving 2 and 2 patients; hyperbilirubinaemia is judged similarly, with three of the four laboratory-rule patients excluded on adjudication (one due to HLH progression, two within assay variability), leaving 1 patient. Hypofibrinogenaemia at 2 weeks (fibrinogen < 2.0 g/L) occurred in 18 patients (75.0%), reported as grade 3-4 in 0: one patient sits at 0.99 g/L, just below the 1.0 g/L CTCAE grade 3 threshold, and the fall was not graded as an adverse event. No pancreatitis, grade 3-4 hyperamylasaemia, immune-related adverse event, permanent treatment discontinuation, or treatment-related death occurred. Asparaginase was omitted in 4 patients (16.7%; protocol-permitted for organ intolerance) and liposomal doxorubicin was given at 25 mg/m² in the 4 patients older than 60 years (16.7%). Every non-haematological row has a grade 3-4 count of 0 because none of those events carries a CTCAE grade in the source; the 17 (70.8%) total is counted from haematological grades alone and is not the sum of the rows. The derivation of every cell is documented in `output/06_table2.csv`.

---

## 6. What this dataset cannot support

The following limitations are stated in the manuscript and apply to any secondary analysis of these data:

1. Single-arm design with no comparator. No attribution of effect to asparaginase, to PD-1 blockade, or to the DEP backbone is possible.
2. Twenty-four patients, eleven of them in histological subgroups of one or two patients. Subgroup estimates are descriptive only.
3. The primary endpoint measures HLH activity, not lymphoma response. Formal lymphoma response was not assessed at a fixed time point.
4. Only 8 patients had a 4-week HLH assessment and 6 a 4-week EBV-DNA measurement, because most responders were bridged to lymphoma-directed therapy.
5. Plasma EBV-DNA assays were not centralised, so absolute values are not comparable between centres; the analyses focus on within-patient change.
6. Cytopenia grades reflect a population in whom cytopenia is partly a manifestation of HLH itself, and transfusion was not corrected for.
7. Median follow-up is 8.5 months, five patients were lost to follow-up, the 12-month survival estimate rests on 9 patients at risk and the 36-month estimate on 2.

---

## 7. Ethical and legal

The study was approved by the Beijing Friendship Hospital Ethics Committee (No. 2023-P2-040-01), conducted in accordance with the Declaration of Helsinki and the study protocol, and all participants gave written informed consent, which included consent for de-identified data to be shared for research purposes.

No re-identification will be attempted, and users of these data must not attempt it.

---

## 7.1 Version history

| Version | Date | Change |
|---|---|---|
| 1.1 | 13 November 2026 | First public release |
| 1.2 | 14 November 2026 | Two changes. (1) `regimen` for pid 1 corrected from `DEP+PD1` to `DEP+PD-1`; no analytical output changed, the correction affects only the literal string used to select patients in whom asparaginase was omitted. (2) `code/06_table2_toxicity.py` no longer fills an event count into the grade 3–4 column of the non-haematological rows. The module docstring had always described those cells as investigator-adjudicated, but the code filled pneumonia and CMV viraemia with the event count (2 and 2) and carried atrial fibrillation as graded 3–4 (1), and counted hypofibrinogenaemia from the laboratory value (1). The source record states that these events occurred and carries no CTCAE grade for them, so all four cells are now 0, matching Table 2 of the manuscript. The "Any grade 3-4 adverse event" total is unchanged at 17 (70.8%), because it is counted from haematological grades alone. The adjudicated values are held in `ADJUDICATED_GRADE_3_4`. No other row of `06_table2.csv` changed. (3) `subtype_WHO_HAEM5` for seven patients renamed to the WHO-HAEM5 entities used in the revised manuscript: patients 2, 11 and 14 are now `PTCL, NOS`, and patients 4, 10 and 15 are now `systemic CAEBVD, T/NK-cell type (category A2-A3)`. The legacy subtype summary table (`05_table3.csv`) is retained for reproducibility but no longer appears in the revised manuscript. No analytical result quoted in the manuscript changed. |

## 8. Licence

- **Data** (`data/`) — Creative Commons Attribution 4.0 International (CC BY 4.0). See `LICENSE-data.txt`.
- **Code** (`code/`) — MIT. See `LICENSE-code.txt`.

---

## 9. How to cite

**Dataset**

> Pi Y, Wang J, Zhou H, Pan X, Wang Z. De-identified patient-level dataset for L-DEP-PD-1 induction therapy in EBV-positive lymphoma-associated hemophagocytic lymphohistiocytosis (v1.2) [dataset]. Zenodo. 2026. doi:10.5281/zenodo.[V1.2-DATASET-DOI]

**Code**

> Pi Y, Wang J, Zhou H, Pan X, Wang Z. LDEP-PD-1_LA-HLH analysis code (v1.2) [computer software]. Zenodo. 2026. doi:10.5281/zenodo.[V1.2-CODE-DOI]

Earlier versions: dataset doi:10.5281/zenodo.22833596 (v1), doi:10.5281/zenodo.24200366 (v1.1); code doi:10.5281/zenodo.22833311 (v1), doi:10.5281/zenodo.24201414 (v1.1). The concept DOI for the code record, which always resolves to the latest version, is doi:10.5281/zenodo.22833310.

---

## 10. Contact

Zhao Wang, MD — Department of Hematology, Beijing Friendship Hospital, Capital Medical University, 95 Yongan Road, Xicheng District, Beijing 100050, China.
E-mail: wangzhao@ccmu.edu.cn

Correspondence about the dataset itself: Yubo Pi (first author).
