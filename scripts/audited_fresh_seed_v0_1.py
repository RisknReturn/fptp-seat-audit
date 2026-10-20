"""Direct check of the audited regional rule at new outer points (review of 2026-10-21, item 3).
Design (written before running): the rule (3,3,3,9,5,3) was chosen with information from the seed-0 points. Here the
rule is fixed and run directly at fresh, independent outer points: 2^12 x 8 scrambled Sobol, seed 3000 (replicate seeds
3000 + 7919 r; not used anywhere before), lambda = 0, the same model and inputs as Table N7. No paired increment from
seed 0 is added. Comparison with the pooled reference of Table N7a exactly as stored in out/r1_note_tables_v02.json
(all 96 scrambles of the three conditional RQMC references, equal weight, SE from the replicate values, df 95),
per-event nominal 95% t approximate intervals with Welch-Satterthwaite degrees of freedom (hybrid: 8 replicates,
df 7), and the same bands as N7 (1e-3, 5e-4). Like N7a, the regional-rule bias is not
subtracted. This is a descriptive confirmation; no pass/fail rule is attached, and the result is reported as it comes.
Writes out/audited_fresh_seed_v0_1.json (with all replicate values)."""
import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import sys, json, time, numpy as np
from scipy.stats import t as student_t
sys.path.insert(0, _ROOT)
src = open(f"{_ROOT}/scripts/v02_experiments.py").read()
sys.argv = [sys.argv[0], "4"]
exec(src[:src.index('out = {"kappa"')])
mp2 = mk(C.POLLS_2024)
EVK = list(events(QEngine2(mp2, nat_on=False, mu_fixed=np.zeros(4)).run(reps=1, m_log2=0)["dist"]).keys())
O = f"{_ROOT}/out"
NT = json.load(open(f"{O}/r1_note_tables_v02.json"))["events"]
t0 = time.time()
r_ = QEngine2(mp2, reg_ns=(3, 3, 3, 9, 5, 3)).run(lam=0.0, m_log2=12, reps=8, seed=3000)
sec = time.time() - t0
E = np.array([[events(x)[k] for k in EVK] for x in r_["dist_reps"]])
h, hs = E.mean(0), E.std(0, ddof=1) / np.sqrt(len(E))
rm = np.array([NT[e]["ref_pooled"] for e in EVK]); rse = np.array([NT[e]["ref_se"] for e in EVK]); rdf = 95.0
d = h - rm; se = np.sqrt(hs ** 2 + rse ** 2)
df = se ** 4 / (hs ** 4 / 7 + rse ** 4 / rdf)
q = student_t.ppf(0.975, df)
lo, hi = d - q * se, d + q * se
rows = {e: dict(hybrid=float(h[i]), hybrid_se=float(hs[i]), ref=float(rm[i]), ref_se=float(rse[i]), diff=float(d[i]),
                ci95=[float(lo[i]), float(hi[i])], df=float(df[i]), within_1e3=bool(lo[i] >= -1e-3 and hi[i] <= 1e-3),
                within_5e4=bool(lo[i] >= -5e-4 and hi[i] <= 5e-4)) for i, e in enumerate(EVK)}
res = dict(design=dict(rule=[3, 3, 3, 9, 5, 3], m_log2=12, reps=8, seed=3000, lam=0.0), sec=sec, reps_events=E.tolist(),
           events=rows, n_within_1e3=int(sum(r["within_1e3"] for r in rows.values())),
           n_within_5e4=int(sum(r["within_5e4"] for r in rows.values())), tv_max_node=r_["tv_max_node"])
json.dump(res, open(f"{O}/audited_fresh_seed_v0_1.json", "w"), indent=1)
print(f"{sec:.0f}s  within 1e-3: {res['n_within_1e3']}/18, within 5e-4: {res['n_within_5e4']}/18", flush=True)
for e, r in rows.items():
    print(f"  {e[:26]:26s} diff {r['diff']:+.1e}  [{r['ci95'][0]:+.1e}, {r['ci95'][1]:+.1e}]  df {r['df']:.0f}", flush=True)
