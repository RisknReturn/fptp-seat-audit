"""R1 review 1.1 / 4.2: lambda-difference precision of the hybrid engine with common random numbers.
Same outer Sobol points (seed 0, 2^11 x 4) at lambda = 0 and 0.5, default regional rule; paired per-replicate
differences. Compared in the R1 response with conditional RQMC's CRN differences (cond_mc_v02.py)."""
import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import sys, json, time, numpy as np
sys.path.insert(0, _ROOT)
src = open(f"{_ROOT}/scripts/v02_experiments.py").read()
sys.argv = [sys.argv[0], "4"]
exec(src[:src.index('out = {"kappa"')])
mp = mk(C.POLLS_2024)
EVK = list(events(QEngine2(mp, nat_on=False, mu_fixed=np.zeros(4)).run(reps=1, m_log2=0)["dist"]).keys())
res = {}
E = {}
for lam in (0.0, 0.5):
    t0 = time.time(); r = QEngine2(mp).run(lam=lam, m_log2=11, reps=4, seed=0); sec = time.time() - t0
    E[lam] = np.array([[events(x)[k] for k in EVK] for x in r["dist_reps"]])
    res[str(lam)] = dict(sec=sec, mean=dict(zip(EVK, E[lam].mean(0).tolist())))
    print(f"lam {lam}: {sec:.0f}s", flush=True)
D = E[0.5] - E[0.0]
res["diff"] = dict(zip(EVK, D.mean(0).tolist()))
res["diff_se_crn"] = dict(zip(EVK, (D.std(0, ddof=1) / np.sqrt(len(D))).tolist()))
json.dump(res, open(f"{_ROOT}/out/hybrid_lam_crn_v02.json", "w"), indent=1)
print("done", {k: (round(res['diff'][k], 5), f"{res['diff_se_crn'][k]:.1e}") for k in EVK}, flush=True)
