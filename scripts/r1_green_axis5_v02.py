"""R1 follow-up 3.5 step 2: axis 5 to 7 points, and the combined rule (3,3,3,9,5,3) to check additivity of the one-axis changes at the same outer
points (2^11 x 4, seed 0), paired with the stored (3,3,3,3,3,3) replicates. Appends to out/r1_followup_v02.json."""
import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import sys, json, time, numpy as np
sys.path.insert(0, _ROOT)
src = open(f"{_ROOT}/scripts/v02_experiments.py").read()
sys.argv = [sys.argv[0], "4"]
exec(src[:src.index('out = {"kappa"')])
mp2 = mk(C.POLLS_2024)
EVK = list(events(QEngine2(mp2, nat_on=False, mu_fixed=np.zeros(4)).run(reps=1, m_log2=0)["dist"]).keys())
rg = json.load(open(f"{_ROOT}/out/reg_convergence_v02b.json"))
E333 = np.array(rg["paired"]["(3, 3, 3, 3, 3, 3)"]["reps"])
OUT = f"{_ROOT}/out/r1_followup_v02.json"
res = json.load(open(OUT))
k = EVK.index("Grn>=4 (actual or more)")
for rule in ((3, 3, 3, 3, 7, 3), (3, 3, 3, 9, 5, 3)):

    t0 = time.time(); r_ = QEngine2(mp2, reg_ns=rule).run(m_log2=11, reps=4, seed=0); sec = time.time() - t0
    E = np.array([[events(x)[kk] for kk in EVK] for x in r_["dist_reps"]])
    Dd = E - E333
    res["green_axes"][str(rule)] = dict(sec=sec, mean=dict(zip(EVK, E.mean(0).tolist())),
                                         diff_vs_333333=dict(zip(EVK, Dd.mean(0).tolist())),
                                         diff_se=dict(zip(EVK, (Dd.std(0, ddof=1) / 2).tolist())))
    print(f"rule {rule}: {sec:.0f}s  dGrn {Dd[:, k].mean():+.3e} (se {Dd[:, k].std(ddof=1)/2:.1e})", flush=True)
    json.dump(res, open(OUT, "w"), indent=1)
print("done", flush=True)
