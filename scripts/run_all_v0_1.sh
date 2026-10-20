#!/bin/bash
cd "$(dirname "$0")/.."   # repository root
set -x
python3 scripts/converge_v0_1.py > out/converge_v0_1.log 2>&1
python3 scripts/mc_compare_v0_1.py > out/mc_compare_std_v0_1.log 2>&1
python3 scripts/backtest_v0_1.py 3,3,2,2,2,2 > out/backtest_v0_1.log 2>&1
python3 scripts/ref_exact_v0_1.py > out/ref_exact_v0_1.log 2>&1
echo ALLDONE
