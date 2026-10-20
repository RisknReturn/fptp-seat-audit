"""mu prior as Student-t predictive from n=3 past transitions (df=2, scale sqrt(4/3)) instead of N(mean, sd)."""
import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import sys, json, numpy as np
sys.path.insert(0, _ROOT)
sys.argv = [sys.argv[0], "4"]
exec(open(f"{_ROOT}/scripts/v02_experiments.py").read().split("out = {")[0])
SC = np.sqrt(4 / 3)
out = {}
mo = mk(POLLS_ACT)
out["oracle_mu_t"] = summ(QEngine2(mo, nat_on=False).run(m_log2=10, reps=4, mu_df=2, mu_scale=SC), "oracle, mu ~ t2 predictive")
mp = mk(C.POLLS_2024)
out["poll_mu_t"] = summ(QEngine2(mp).run(m_log2=11, reps=4, mu_df=2, mu_scale=SC), "final polls, lam=0, mu ~ t2 predictive")
for k in out: out[k].pop("dist", None)
json.dump(out, open(f"{ROOT}/out/mu_tpred_v02.json", "w"), indent=1)
