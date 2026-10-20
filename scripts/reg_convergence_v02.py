"""R1 review 4.1: regional-rule convergence at FIXED outer RQMC points (paired), all 18 events, plus outer-only
expansion at a fixed regional rule. Separates (a) outer scramble SE, (b) regional-rule change, (c) truncation TV.
Model: v0.2.1, final polls, lambda = 0, mu ~ N(2010-19 prior)."""
import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import sys, json, time, numpy as np
sys.path.insert(0, _ROOT)
src = open(f"{_ROOT}/scripts/v02_experiments.py").read()
sys.argv = [sys.argv[0], "4"]
exec(src[:src.index('out = {"kappa"')])
mp = mk(C.POLLS_2024)
RULES = [(3, 3, 2, 2, 2, 2), (3, 3, 3, 3, 3, 3), (5, 5, 2, 2, 2, 2), (5, 5, 3, 3, 3, 3), (5, 5, 5, 3, 3, 3), (5, 5, 5, 5, 3, 3)]
OUTER = [dict(m_log2=11, reps=4, seed=0)]
res = {"rules": {}, "outer": {}}
def ev_reps(r):
    return np.array([[events(x)[k] for k in EVK] for x in r["dist_reps"]])
EVK = list(events(QEngine2(mp, nat_on=False, mu_fixed=np.zeros(4)).run(reps=1, m_log2=0)["dist"]).keys())
base = None
for rule in RULES:
    t0 = time.time(); r = QEngine2(mp, reg_ns=rule).run(**OUTER[0]); sec = time.time() - t0
    E = ev_reps(r); mean = E.mean(0)
    row = dict(rule=list(rule), sec=sec, tv_max_node=r["tv_max_node"], mean=dict(zip(EVK, mean.tolist())),
               se=dict(zip(EVK, (E.std(0, ddof=1) / np.sqrt(len(E))).tolist())))
    if base is None:
        base = E
    else:
        D = E - base   # paired difference against the default rule, same outer points
        row["diff_vs_default"] = dict(zip(EVK, D.mean(0).tolist())); row["diff_se"] = dict(zip(EVK, (D.std(0, ddof=1) / np.sqrt(len(D))).tolist()))
    res["rules"][str(rule)] = row
    big = {k: round(v, 6) for k, v in row.get("diff_vs_default", {}).items() if abs(v) > 4e-4}
    print(f"rule {rule}: {sec:.0f}s  Grn {row['mean']['Grn>=4 (actual or more)']:.6f} PC {row['mean']['PC>=4 (actual or more)']:.6f} "
          f"Oth {row['mean']['Oth>=5 (actual or more)']:.6f}  |diff|>4e-4: {big}", flush=True)
    json.dump(res, open(f"{_ROOT}/out/reg_convergence_v02.json", "w"), indent=1)
for cfg in (dict(m_log2=12, reps=8, seed=0), dict(m_log2=12, reps=16, seed=100), dict(m_log2=13, reps=16, seed=200)):
    t0 = time.time(); r = QEngine2(mp).run(**cfg); sec = time.time() - t0
    E = ev_reps(r)
    row = dict(cfg=cfg, sec=sec, mean=dict(zip(EVK, E.mean(0).tolist())), se=dict(zip(EVK, (E.std(0, ddof=1) / np.sqrt(len(E))).tolist())),
               oth_reps=E[:, EVK.index("Oth>=5 (actual or more)")].tolist())
    res["outer"][str(cfg)] = row
    print(f"outer {cfg}: {sec:.0f}s  Oth {row['mean']['Oth>=5 (actual or more)']:.6f} (se {row['se']['Oth>=5 (actual or more)']:.6f}) "
          f"Grn {row['mean']['Grn>=4 (actual or more)']:.6f}", flush=True)
    json.dump(res, open(f"{_ROOT}/out/reg_convergence_v02.json", "w"), indent=1)
print("done")
