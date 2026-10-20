"""R1 review 4.1, part 2 (resume after the first run hit the time limit during rule (5,5,5,5,3,3)).
Keeps the five finished rules from out/reg_convergence_v02.json, then
  (a) re-runs the default rule and (3,3,3,3,3,3) at the same fixed outer points to get paired per-replicate arrays,
  (b) adds (3,3,3,5,5,5): sensitivity of regional dims 4-6 (the first run varied dims 1-3 only),
  (c) outer-only expansions at the default regional rule.
Model: v0.2.1, final polls, lambda = 0, mu ~ N(2010-19 prior)."""
import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import sys, json, time, numpy as np
sys.path.insert(0, _ROOT)
src = open(f"{_ROOT}/scripts/v02_experiments.py").read()
sys.argv = [sys.argv[0], "4"]
exec(src[:src.index('out = {"kappa"')])
mp = mk(C.POLLS_2024)
OUT = f"{_ROOT}/out/reg_convergence_v02b.json"
res = json.load(open(f"{_ROOT}/out/reg_convergence_v02.json"))
res["note"] = "rules from run 1 kept; rule (5,5,5,5,3,3) of run 1 not finished (time limit); paired arrays re-run here"
res["paired"] = {}; res["outer"] = {}
OUTER = dict(m_log2=11, reps=4, seed=0)
EVK = list(events(QEngine2(mp, nat_on=False, mu_fixed=np.zeros(4)).run(reps=1, m_log2=0)["dist"]).keys())
def ev_reps(r):
    return np.array([[events(x)[k] for k in EVK] for x in r["dist_reps"]])
E = {}
for rule in [(3, 3, 2, 2, 2, 2), (3, 3, 3, 3, 3, 3), (3, 3, 3, 5, 5, 5)]:
    t0 = time.time(); r = QEngine2(mp, reg_ns=rule).run(**OUTER); sec = time.time() - t0
    E[rule] = ev_reps(r)
    res["paired"][str(rule)] = dict(sec=sec, tv_max_node=r["tv_max_node"], reps=E[rule].tolist(),
                                    mean=dict(zip(EVK, E[rule].mean(0).tolist())))
    print(f"rule {rule}: {sec:.0f}s  Grn {E[rule][:, EVK.index('Grn>=4 (actual or more)')].mean():.6f}", flush=True)
    json.dump(res, open(OUT, "w"), indent=1)
def pdiff(a, b):
    D = E[a] - E[b]
    return dict(mean=dict(zip(EVK, D.mean(0).tolist())), se=dict(zip(EVK, (D.std(0, ddof=1) / np.sqrt(len(D))).tolist())))
res["paired_diffs"] = {"(3,3,3,3,3,3)-(3,3,2,2,2,2)": pdiff((3, 3, 3, 3, 3, 3), (3, 3, 2, 2, 2, 2)),
                       "(3,3,3,5,5,5)-(3,3,3,3,3,3)": pdiff((3, 3, 3, 5, 5, 5), (3, 3, 3, 3, 3, 3))}
d = res["paired_diffs"]["(3,3,3,5,5,5)-(3,3,3,3,3,3)"]
print("dims 4-6 sensitivity, |diff|>1e-4:", {k: round(v, 6) for k, v in d["mean"].items() if abs(v) > 1e-4}, flush=True)
json.dump(res, open(OUT, "w"), indent=1)
for cfg in (dict(m_log2=12, reps=8, seed=0), dict(m_log2=12, reps=16, seed=100), dict(m_log2=13, reps=16, seed=200)):
    t0 = time.time(); r = QEngine2(mp).run(**cfg); sec = time.time() - t0
    Eo = ev_reps(r)
    row = dict(cfg=cfg, sec=sec, mean=dict(zip(EVK, Eo.mean(0).tolist())), se=dict(zip(EVK, (Eo.std(0, ddof=1) / np.sqrt(len(Eo))).tolist())),
               oth_reps=Eo[:, EVK.index("Oth>=5 (actual or more)")].tolist())
    res["outer"][str(cfg)] = row
    print(f"outer {cfg}: {sec:.0f}s  Oth {row['mean']['Oth>=5 (actual or more)']:.6f} (se {row['se']['Oth>=5 (actual or more)']:.6f}) "
          f"Grn {row['mean']['Grn>=4 (actual or more)']:.6f}", flush=True)
    json.dump(res, open(OUT, "w"), indent=1)
print("done")
