"""Poisson-binomial kernel benchmark only (review of 2026-10-21, issue 2): dense recursion vs the active-support kernel
with the hybrid's near-certain truncation (eps = 1e-12), named EXPLICITLY (pb_dense, pb_active from cond_mc_v02.py),
so the comparison does not depend on UKSEAT_PB. Same design as part (1) of r1_followup_v02.py: the same 1,024 factor
points (scrambled Sobol, seed default_rng([20261030, 0, 1])), model and full PMFs, one lambda and two lambdas (common
random numbers), warm runs (JIT outside timing), median of 3. Records each kernel's function name.
Writes out/pb_kernel_bench_v0_1.json; does not touch out/r1_followup_v02.json (whose part (2) is a separate audit).
Timings depend on the machine and on other load; the speed ratio, not the seconds, is the quantity of interest."""
import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import sys, json, time, platform, numpy as np
import numba
from numba import njit, prange
from scipy.stats import qmc, norm
sys.path.insert(0, _ROOT)
src = open(f"{_ROOT}/scripts/cond_mc_v02.py").read()
sys.argv = [sys.argv[0], "rqmc", "1024", "1", "0"]
exec(src[:src.index("EVK = None")].replace("cache=True", "cache=False"))   # model, pb_dense, pb_active, one_rep, shift
assert pb_dense is not pb_active and pb_dense.__name__ == "pb_dense" and pb_active.__name__ == "pb_active"

Z = norm.ppf(np.clip(qmc.Sobol(D, scramble=True, seed=np.random.default_rng([20261030, 0, 1])).random_base2(10), 1e-12, 1 - 1e-12))
all_shift = dict(shift)
bench, results = [], {}
for name, kern in (("dense", pb_dense), ("active_eps1e-12", pb_active)):
    globals()["pb_all"] = kern                  # one_rep() looks up pb_all at call time
    for ls in ((0.0,), (0.0, 0.5)):
        shift = {x: all_shift[x] for x in ls}
        one_rep(Z[:4])                          # warm-up / JIT outside timing
        ts = []
        for _ in range(3):
            t0 = time.perf_counter(); r_ = one_rep(Z); ts.append(time.perf_counter() - t0)
        results[(name, ls)] = r_
        bench.append(dict(kernel=name, function=kern.__name__, lambdas=list(ls), N=len(Z), seconds=ts,
                          median=float(np.median(ts))))
        print("bench", bench[-1], flush=True)
shift = all_shift
res = dict(pb_benchmark=bench,
           pb_max_pmf_diff=float(max(np.abs(results[("dense", ls)][lam] - results[("active_eps1e-12", ls)][lam]).max()
                                     for ls in ((0.0,), (0.0, 0.5)) for lam in ls)),
           pb_tv_bound_per_party=n * 1e-12,
           speed_ratio_dense_over_active={str(ls): bench[i]["median"] / bench[i + 2]["median"]
                                         for i, ls in enumerate(((0.0,), (0.0, 0.5)))},
           environment=dict(python=platform.python_version(), numpy=np.__version__, numba=numba.__version__,
                            numba_threads=numba.get_num_threads(), machine=platform.machine(), system=platform.system()))
st = json.load(open(f"{_ROOT}/out/r1_followup_v02.json"))
res["stored_r1_followup_v02"] = dict(pb_max_pmf_diff=st["pb_max_pmf_diff"],
                                     speed_ratio_dense_over_active={str(tuple(st["pb_benchmark"][i]["lambdas"])):
                                         st["pb_benchmark"][i]["median"] / st["pb_benchmark"][i + 2]["median"]
                                         for i in range(2)})
json.dump(res, open(f"{_ROOT}/out/pb_kernel_bench_v0_1.json", "w"), indent=1)
print("max PMF diff", res["pb_max_pmf_diff"], " ratios", res["speed_ratio_dense_over_active"],
      " stored", res["stored_r1_followup_v02"], flush=True)
