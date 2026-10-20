#!/bin/bash
cd "$(dirname "$0")/.."   # repository root
python3 scripts/ref_exact2_v0_1.py > out/ref_exact2_v0_1.log 2>&1
python3 scripts/mc_compare_v0_1.py > out/mc_compare_std_v0_1.log 2>&1
python3 scripts/backtest2_v0_1.py 3,3,2,2,2,2 > out/backtest2_v0_1.log 2>&1
python3 scripts/tables_v0_1.py > out/tables_v0_1.log 2>&1
python3 scripts/plots_v0_1.py > out/plots_v0_1.log 2>&1
echo ALLDONE2
