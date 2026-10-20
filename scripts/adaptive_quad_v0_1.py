"""Goal-oriented (dimension-adaptive) choice of quadrature points for the seat engine.

Greedy anisotropic refinement: at each step, refine by one level the dimension whose refinement changes the
target quantities (18 event probabilities) the most per unit of compute time; stop when no single refinement
changes any target by more than TOL. Dimensions: 5 national factors (party-normalised impact rotation; dim 0 uses
the trapezoid rule) and 6 regional factor indices (same rule in every region). Validated against plain MC (3e6).
"""
import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import sys, json, time, numpy as np
sys.path.insert(0, _ROOT)
from ukseat.setup2024 import make_model
from ukseat.engine import FastEngine
from ukseat.report import events

TOL = float(sys.argv[1]) if len(sys.argv) > 1 else 4e-4
BUDGET = float(sys.argv[2]) if len(sys.argv) > 2 else 3600.0
NAT0 = [("trap", 0.7, 7.5), ("trap", 0.5, 7.5), ("trap", 0.35, 7.5), ("trap", 0.25, 7.5)]
NATK = [("gh", 2), ("gh", 3), ("gh", 5), ("gh", 7), ("gh", 9)]
REGL = [2, 3, 5, 7]
m = make_model()
mc = json.load(open(f"{_ROOT}/out/mc_compare_std_v0_1.json"))
m6 = [r for r in mc["runs"] if r["N"] == 1_000_000]
mcev = {k: float(np.mean([r["ev"][k] for r in m6])) for k in m6[0]["ev"]}
cache = {}

def spec(lv):
    nat = [NAT0[lv[0]]] + [NATK[l] for l in lv[1:5]]
    reg = tuple(REGL[l] for l in lv[5:11])
    return nat, reg

def run(lv):
    lv = tuple(lv)
    if lv in cache:
        return cache[lv]
    nat, reg = spec(lv)
    r = FastEngine(m, nat_spec=nat, reg_ns=reg, nat_norm=True).run(0.0)
    ev = events(r["dist"])
    cache[lv] = dict(ev=ev, sec=r["seconds"], nodes=r["nodes_outer"], reg=r["nodes_reg"])
    return cache[lv]

def gap(a, b):
    return max(abs(a[k] - b[k]) for k in a)

maxlv = [len(NAT0) - 1] + [len(NATK) - 1] * 4 + [len(REGL) - 1] * 6
lv = [0] * 11
run(lv)  # compile + base
log = []
t0 = time.time()
while time.time() - t0 < BUDGET:
    cur = run(lv)
    cands = []
    for d in range(11):
        if lv[d] >= maxlv[d]:
            continue
        nl = list(lv); nl[d] += 1
        r = run(nl)
        eta = gap(r["ev"], cur["ev"])
        cands.append((eta / max(r["sec"], 1e-3), eta, d, nl, r))
    if not cands:
        break
    best = max(cands, key=lambda c: c[0])
    etamax = max(c[1] for c in cands)
    errmc = gap(cur["ev"], mcev)
    step = dict(levels=list(lv), spec=str(spec(lv)), nodes=cur["nodes"], sec=cur["sec"], err_vs_mc=errmc,
                eta={int(c[2]): c[1] for c in cands}, chosen=best[2], eta_max=etamax, elapsed=time.time() - t0)
    log.append(step)
    print(f"[{time.time()-t0:6.0f}s] levels={lv} nodes={cur['nodes']:6d} sec={cur['sec']:6.1f} "
          f"max|ev-MC|={errmc:.4f}  eta_max={etamax:.1e} -> refine dim {best[2]} (eta {best[1]:.1e})", flush=True)
    json.dump(dict(log=log, mc=mcev, tol=TOL), open(f"{_ROOT}/out/adaptive_quad_v0_1.json", "w"))
    if etamax < TOL:
        print("converged: no single refinement changes any event by more than", TOL, flush=True)
        break
    lv = best[3]
final = run(lv)
print("final", spec(lv), "nodes", final["nodes"], "sec %.1f" % final["sec"], "max|ev-MC| %.4f" % gap(final["ev"], mcev))
worst = sorted(((abs(final["ev"][k] - mcev[k]), k) for k in mcev), reverse=True)[:5]
print("largest gaps vs MC:", [(k, round(v, 4)) for v, k in worst])
json.dump(dict(log=log, mc=mcev, tol=TOL, final=dict(levels=lv, spec=str(spec(lv)), **{k: final[k] for k in ("sec", "nodes", "reg")},
               ev=final["ev"])), open(f"{_ROOT}/out/adaptive_quad_v0_1.json", "w"))
