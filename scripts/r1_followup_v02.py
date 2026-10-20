"""R1 follow-up review: (1) conditional PB kernel with the same near-certain truncation as the hybrid (eps = 1e-12)
vs the dense kernel, same factor points / model / full PMFs, one lambda and two lambdas (CRN), warm runs, median of 3;
(2) Green: which of regional axes 4-6 drives the +6.89e-4 change of (3,3,3,5,5,5) vs (3,3,3,3,3,3): one axis at a time,
same outer points (2^11 x 4, seed 0), paired with the stored (3,3,3,3,3,3) replicates.
Writes out/r1_followup_v02.json.
Review of 2026-10-21: the dense side of (1) was `pb_all`, which follows UKSEAT_PB and has defaulted to pb_active since
v0.2.1, so a default re-run compared active with active. It is now pb_dense explicitly and each entry records the
function name. The stored out/r1_followup_v02.json was produced while the default was still dense (its max PMF
difference 1.4e-14 is between two different kernels). To re-check (1) alone, run scripts/pb_kernel_bench_v0_1.py;
re-running this file also re-runs (2) and overwrites the stored json."""
import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import sys, json, time, numpy as np
from numba import njit, prange
from scipy.stats import qmc, norm
sys.path.insert(0, _ROOT)
src = open(f"{_ROOT}/scripts/cond_mc_v02.py").read()
sys.argv = [sys.argv[0], "rqmc", "1024", "1", "0"]
exec(src[:src.index("EVK = None")].replace("cache=True", "cache=False"))   # model, pb_all (dense), one_rep, shift
OUT = f"{_ROOT}/out/r1_followup_v02.json"
res = {}


@njit(parallel=True, cache=True)
def pb_active(P):
    n_, K_ = P.shape
    out_ = np.zeros((K_, n_ + 1))
    for k in prange(K_):
        act = np.empty(n_); na = 0; fixed = 0
        for i in range(n_):
            q = P[i, k]
            if q <= 1e-12:
                continue
            if q >= 1.0 - 1e-12:
                fixed += 1
            else:
                act[na] = q; na += 1
        d = np.zeros(na + 1); d[0] = 1.0
        for j in range(na):
            q = act[j]
            for s_ in range(j + 1, 0, -1):
                d[s_] = (1 - q) * d[s_] + q * d[s_ - 1]
            d[0] *= 1 - q
        for s_ in range(na + 1):
            out_[k, s_ + fixed] = d[s_]
    return out_


Z = norm.ppf(np.clip(qmc.Sobol(D, scramble=True, seed=np.random.default_rng([20261030, 0, 1])).random_base2(10), 1e-12, 1 - 1e-12))
all_shift = dict(shift)
bench, results = [], {}
dense = pb_dense   # explicit: pb_all follows UKSEAT_PB (default pb_active since v0.2.1)
for name, kern in (("dense", dense), ("active_eps1e-12", pb_active)):
    globals()["pb_all"] = kern
    for ls in ((0.0,), (0.0, 0.5)):
        shift = {x: all_shift[x] for x in ls}
        one_rep(Z[:4])                          # warm-up / JIT outside timing
        ts = []
        for _ in range(3):
            t0 = time.perf_counter(); r_ = one_rep(Z); ts.append(time.perf_counter() - t0)
        results[(name, ls)] = r_
        bench.append(dict(kernel=name, function=kern.__name__, lambdas=list(ls), N=len(Z), seconds=ts, median=float(np.median(ts))))
        print("bench", bench[-1], flush=True)
globals()["pb_all"] = dense
shift = all_shift
res["pb_benchmark"] = bench
res["pb_max_pmf_diff"] = float(max(np.abs(results[("dense", ls)][lam] - results[("active_eps1e-12", ls)][lam]).max()
                                   for ls in ((0.0,), (0.0, 0.5)) for lam in ls))
res["pb_tv_bound_per_party"] = n * 1e-12
json.dump(res, open(OUT, "w"), indent=1)

# (2) one-axis Green audit
src2 = open(f"{_ROOT}/scripts/v02_experiments.py").read()
sys.argv = [sys.argv[0], "4"]
exec(src2[:src2.index('out = {"kappa"')])
mp2 = mk(C.POLLS_2024)
EVK = list(events(QEngine2(mp2, nat_on=False, mu_fixed=np.zeros(4)).run(reps=1, m_log2=0)["dist"]).keys())
rg = json.load(open(f"{_ROOT}/out/reg_convergence_v02b.json"))
E333 = np.array(rg["paired"]["(3, 3, 3, 3, 3, 3)"]["reps"])
res["green_axes"] = {}
for rule in ((3, 3, 3, 5, 3, 3), (3, 3, 3, 3, 5, 3), (3, 3, 3, 3, 3, 5)):
    t0 = time.time(); r_ = QEngine2(mp2, reg_ns=rule).run(m_log2=11, reps=4, seed=0); sec = time.time() - t0
    E = np.array([[events(x)[k] for k in EVK] for x in r_["dist_reps"]])
    Dd = E - E333
    res["green_axes"][str(rule)] = dict(sec=sec, mean=dict(zip(EVK, E.mean(0).tolist())),
                                         diff_vs_333333=dict(zip(EVK, Dd.mean(0).tolist())),
                                         diff_se=dict(zip(EVK, (Dd.std(0, ddof=1) / 2).tolist())))
    k = EVK.index("Grn>=4 (actual or more)"); kp = EVK.index("PC>=4 (actual or more)")
    print(f"rule {rule}: {sec:.0f}s  dGrn {Dd[:, k].mean():+.2e} (se {Dd[:, k].std(ddof=1)/2:.1e})  dPC {Dd[:, kp].mean():+.2e}", flush=True)
    json.dump(res, open(OUT, "w"), indent=1)
print("done", flush=True)
