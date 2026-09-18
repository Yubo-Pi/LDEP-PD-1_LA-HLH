# L-DEP-PD-1 LA-HLH analysis code

Analysis code for the EBV-positive LA-HLH study NCT05775705. Patient-level
data and reference outputs are distributed separately in the companion
Zenodo data package (24 subjects, 118 variables).

## Contents

- `code/02_analysis.py`: statistical analysis.
- `code/run_all.py`: analysis runner; optional private workbook preparation.
- `code/01_data_preparation.py`: source column mapping and preparation.
- `requirements.txt`: Python dependencies with fixed versions.
- `requirements-preparation.txt`: optional workbook-reading dependency.
- `data/README.md`: dataset placement instructions.
- `LICENSE`: CC BY 4.0 license terms.
- `SHA256SUMS.txt`: checksums of release files.

## Reproduce the results

1. Use Python 3.12. Create and activate a virtual environment.
2. Download the companion Zenodo data package and copy its `data/` folder
   into this repository. The required path is
   `data/LDEP_PD1_LA-HLH_deidentified.csv`.
3. From the repository root, run:

```bash
python -m pip install -r requirements.txt
python code/run_all.py
```

`python code/02_analysis.py` is an equivalent analysis entry point once the
CSV has been placed correctly. The analysis creates six result CSVs and
`outputs/reproducibility_check.txt`. It does not regenerate the preparation
log or the data dictionary. Compare results with the Zenodo `outputs/` folder.
The full console log and `table3_patient_level.csv` include subject-level
information; generated data and outputs are excluded by `.gitignore`.

## Optional: prepare from the private workbook

```bash
python -m pip install -r requirements-preparation.txt
python code/run_all.py /path/to/source.xls
```

This requires the authors' private workbook, first sheet, **116 columns**,
two header rows, with the subsequent-therapy free-text column already removed.
The script rejects other column counts and checks selected header anchors.
The source workbook is not included in the public release.

## Companion data and citation

The Zenodo DOI/URL and GitHub repository URL have not yet been assigned.
Add the actual record links to both READMEs after creating the records.
Please cite the associated manuscript and the companion dataset record.

## Associated manuscript

Pi Y, Wang J, Zhou H, Pan X, Wang Z. *Adding asparaginase and PD-1 blockade to the
DEP regimen achieves rapid disease control and Epstein–Barr virus clearance in
EBV-positive lymphoma-associated hemophagocytic lymphohistiocytosis.*
ClinicalTrials.gov identifier: NCT05775705.

## Study design

Prospective, multicentre, single-arm Phase 2 study conducted at three centres in
China. Twenty-four adults with EBV-positive LA-HLH (EBV-positive lymphoma or grade
2–3 EBV-positive T/NK-cell lymphoproliferative disease, meeting HLH-2004 criteria)
received the L-DEP-PD-1 regimen as induction therapy. Patients were enrolled
between August 2023 and July 2026.

Treatment (one cycle, repeated every 2 weeks for up to two cycles):
liposomal doxorubicin 35 mg/m² and etoposide 75 mg/m² on day 1; methylprednisolone
1.5 mg/kg days 1–3 tapered to 0.25 mg/kg by day 14; asparaginase 6,000 IU/m² on
days 2 and 4; camrelizumab 200 mg on day 5. Four patients had asparaginase omitted
because of organ intolerance and are recorded as `regimen = DEP+PD-1`.

The primary endpoint was the overall response rate (ORR: complete or partial
response) two weeks after treatment, graded using the quantifiable HLH markers
proposed by Marsh et al.

## Statistical conventions

The analysis implements the following conventions:

* continuous variables — medians with ranges;
* paired baseline-to-week-2 comparisons — two-sided **asymptotic** Wilcoxon
  signed-rank test;
* paired plasma EBV-DNA — zero-difference pairs excluded, **exact** Wilcoxon
  signed-rank test (n = 16);
* proportions — **Clopper-Pearson exact** 95% confidence intervals;
* overall survival — Kaplan-Meier; median follow-up as the median of individual
  follow-up times;
* two-sided P < 0.05 considered significant.

## Contact

Zhao Wang — Department of Hematology, Beijing Friendship Hospital, Capital
Medical University — wangzhao@ccmu.edu.cn
