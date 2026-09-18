#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_all.py
----------
Runs the complete reproducibility pipeline:

    1. 01_data_preparation.py   (optional; needs the source .xls workbook)
    2. 02_analysis.py           (reproduces all reported results)

Usage
-----
    python code/run_all.py                    # analysis only (uses CSV downloaded from Zenodo)
    python code/run_all.py /path/to/raw.xls   # rebuild CSV, then analyse

Download the companion Zenodo data package and copy its data/ folder into
this repository before running. The source workbook is not publicly included.
"""
from __future__ import annotations

import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable

raw = sys.argv[1] if len(sys.argv) > 1 else None


def run(script, *args):
    cmd = [PY, os.path.join(HERE, script)] + list(args)
    print('\n$ %s' % ' '.join(cmd))
    r = subprocess.run(cmd)
    if r.returncode != 0:
        sys.exit('step failed: %s' % script)


if raw:
    if not os.path.exists(raw):
        sys.exit('raw workbook not found: %s' % raw)
    run('01_data_preparation.py', raw)
else:
    csv = os.path.join(os.path.dirname(HERE), 'data',
                       'LDEP_PD1_LA-HLH_deidentified.csv')
    if not os.path.exists(csv):
        sys.exit('dataset not found: %s\n'
                 'Copy data/LDEP_PD1_LA-HLH_deidentified.csv from the companion Zenodo package into this repository.' % csv)
    print('using downloaded dataset: %s' % os.path.basename(csv))
    print('(pass the path to the source .xls to rebuild it from scratch)')

run('02_analysis.py')
print('\nDone. See outputs/ for all result files.')
