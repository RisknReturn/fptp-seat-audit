"""Where the remaining oracle gap sits: expected seats by nation (England/Scotland/Wales) vs actual, v0.2 realised mu."""
import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import sys, json, numpy as np, pandas as pd
sys.path.insert(0, _ROOT)
sys.argv = [sys.argv[0], "4"]
src = open(f"{_ROOT}/scripts/v02_experiments.py").read()
exec(src[:src.index("out = {\"kappa\"")])
mo = mk(POLLS_ACT)
with Raked(mo, MU_REAL):
    P = seat_marginals(mo, 0.0, False, MU_REAL, ndraw=600)
P0 = seat_marginals(mo, 0.0, False, np.zeros(4), ndraw=600)
lab = ["Con", "Lab", "LD", "Ref", "Grn", "SNP", "PC", "OthMax"]
act = mo.df[[f"a24_{k}" for k in lab]].values.argmax(1)
A = np.zeros_like(P); A[np.arange(len(act)), act] = 1
cty = mo.df.country.values
rows = []
for c in ["England", "Scotland", "Wales"]:
    i = cty == c
    for nm, M in (("v0.2 mu=0", P0), ("v0.2 mu=real", P), ("actual", A)):
        rows.append(dict(nation=c, run=nm, **{k: round(float(M[i, j].sum()), 1) for j, k in enumerate(CATS)}))
t = pd.DataFrame(rows); print(t.to_string(index=False))
t.to_csv(f"{ROOT}/out/oracle_by_country_v02.csv", index=False)
