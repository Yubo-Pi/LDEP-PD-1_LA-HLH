"""
run_all.py — execute the full analysis in order.

    cd LDEP-PD-1_LA-HLH_v1.1
    python code/run_all.py

Every script writes human-readable text and CSV output to output/.
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = [
    '01_baseline_and_response.py',
    '02_ebv.py',
    '03_laboratory_changes.py',
    '04_survival.py',
    '05_tables_3_4.py',
    '06_table2_toxicity.py',
]

env = dict(os.environ)
env['PYTHONPATH'] = HERE + os.pathsep + env.get('PYTHONPATH', '')

failed = []
for s in SCRIPTS:
    print('\n' + '#' * 78)
    print('# ' + s)
    print('#' * 78)
    r = subprocess.run([sys.executable, os.path.join(HERE, s)], env=env)
    if r.returncode != 0:
        failed.append(s)

print('\n' + '=' * 78)
if failed:
    print('FAILED:', ', '.join(failed))
    sys.exit(1)
print('All scripts completed. See output/ for results.')
