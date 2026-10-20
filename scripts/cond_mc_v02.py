"""R1 review 4.2: conditional MC and conditional RQMC baselines.

All common factors F = (z: 5, mu: 4, regional w: 57) are sampled (66 dims); given F, seats are independent and each
party's seat PMF is computed exactly by the Poisson-binomial recursion (Rao-Blackwellised: no seat-level sampling).
Same model as the hybrid engine (v0.2.1, final polls, mu ~ N(2010-19 prior)); no regional Gauss-Hermite, so no
regional-rule bias. Common random numbers: lambda = 0 and lambda = 0.5 use the same F draws.
usage: python cond_mc_v02.py MODE N REPS SEED   (MODE = iid | rqmc; N = samples per replicate, power of 2 for rqmc)
PB kernel: active-support with eps = 1e-12 by default (v0.2.1); UKSEAT_PB=dense reproduces the stored v0.2 runs.
Seeds: replicate r uses default_rng([SEED, r, 1]) (Sobol scramble) or default_rng([SEED, r, 2]) (iid), so the
draws do not overlap with the reviewer's audit (which used SEED + 1009 r).
"""
import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import sys, json, time, numpy as np
from numba import njit, prange
from scipy.stats import qmc, norm
sys.path.insert(0, _ROOT)
MODE, N, REPS, SEED = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
src = open(f"{_ROOT}/scripts/v02_experiments.py").read()
sys.argv = [sys.argv[0], "4"]
exec(src[:src.index('out = {"kappa"')])
mp = mk(C.POLLS_2024)
K = len(CATS); n = mp.n
pres = np.isfinite(mp.logb[:, :K])
base0 = np.where(pres, mp.logb[:, :K], 0.0) + mp.mean_delta()[:, :K]
shift = {0.0: mp.scenario_shift(0.0)[:K], 0.5: mp.scenario_shift(0.5)[:K]}
A = mp.A[:K]; L = mp.L[:, :K, :]
regs = [(np.where(mp.reg == j)[0], Cj[:K]) for j, Cj in enumerate(mp.C)]
dz, dm = A.shape[1], L.shape[2]; dr = sum(Cj.shape[1] for _, Cj in regs); D = dz + dm + dr
print(f"dims: national {dz} + mu {dm} + regional {dr} = {D}; seats {n}; regions {len(regs)}", flush=True)

@njit(parallel=True, cache=True)
def pb_dense(P):
    """v0.2 kernel used for the stored reference runs: every q > 0 enters the dense recursion."""
    n, K = P.shape
    out = np.zeros((K, n + 1))
    for k in prange(K):
        d = np.zeros(n + 1); d[0] = 1.0; m = 0
        for i in range(n):
            q = P[i, k]
            if q <= 0.0:
                continue
            m += 1
            for s in range(m, 0, -1):
                d[s] = d[s] * (1 - q) + d[s - 1] * q
            d[0] *= 1 - q
        out[k] = d
    return out


@njit(parallel=True, cache=True)
def pb_active(P):
    """Same near-certain truncation as the hybrid (eps = 1e-12): q <= eps dropped, q >= 1 - eps counted as fixed seats.
    Per-party TV change <= n * eps (631e-12). R1 follow-up: about 2.3-2.4x faster than pb_dense; max PMF difference
    1.5e-14 on 1,024 factor points."""
    n, K = P.shape
    out = np.zeros((K, n + 1))
    for k in prange(K):
        act = np.empty(n); na = 0; fixed = 0
        for i in range(n):
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
            for s in range(j + 1, 0, -1):
                d[s] = (1 - q) * d[s] + q * d[s - 1]
            d[0] *= 1 - q
        for s in range(na + 1):
            out[k, s + fixed] = d[s]
    return out


import os as _os2
pb_all = pb_dense if _os2.environ.get("UKSEAT_PB", "active") == "dense" else pb_active   # v0.2.1 default: active

def probs(eta):
    return mp.win_probs(eta, pres)

def one_rep(Zs):
    acc = {lam: np.zeros((K, n + 1)) for lam in shift}
    for Zrow in Zs:
        z = Zrow[:dz]; mu = mp.mu_mean + mp.mu_sd * Zrow[dz:dz + dm]; off = dz + dm
        eta = base0 + (A @ z)[None, :] + np.einsum("ikj,j->ik", L, mu)
        for idx, Cj in regs:
            w = Zrow[off: off + Cj.shape[1]]; off += Cj.shape[1]
            eta[idx] += (Cj @ w)[None, :]
        for lam, sh in shift.items():
            acc[lam] += pb_all(probs(eta + sh[None, :]))
    return {lam: a / len(Zs) for lam, a in acc.items()}

EVK = None; out = dict(mode=MODE, N=N, reps=REPS, seed=SEED, dims=D, reps_events={}, seconds=[],
                      pb_kernel=_os2.environ.get("UKSEAT_PB", "active"))
t_all = time.time()
for r in range(REPS):
    t0 = time.time()
    if MODE == "rqmc":
        U = qmc.Sobol(D, scramble=True, seed=np.random.default_rng([SEED, r, 1])).random(N)
        Zs = norm.ppf(np.clip(U, 1e-12, 1 - 1e-12))
    else:
        Zs = np.random.default_rng([SEED, r, 2]).standard_normal((N, D))
    dist = one_rep(Zs)
    sec = time.time() - t0; out["seconds"].append(sec)
    for lam, d in dist.items():
        ev = events(d); EVK = list(ev)
        out["reps_events"].setdefault(str(lam), []).append([ev[k] for k in EVK])
    e0 = out["reps_events"]["0.0"][-1]
    print(f"rep {r}: {sec:.1f}s  Grn {e0[EVK.index('Grn>=4 (actual or more)')]:.6f} PC {e0[EVK.index('PC>=4 (actual or more)')]:.6f} "
          f"Oth {e0[EVK.index('Oth>=5 (actual or more)')]:.6f} Lab326 {e0[EVK.index('Lab>=326 (majority)')]:.6f}", flush=True)
out["events"] = EVK
for lam in out["reps_events"]:
    E = np.array(out["reps_events"][lam])
    out[f"mean_{lam}"] = dict(zip(EVK, E.mean(0).tolist())); out[f"se_{lam}"] = dict(zip(EVK, (E.std(0, ddof=1) / np.sqrt(len(E))).tolist()))
Dd = np.array(out["reps_events"]["0.5"]) - np.array(out["reps_events"]["0.0"])
out["diff_lam05_minus_0"] = dict(zip(EVK, Dd.mean(0).tolist())); out["diff_se_crn"] = dict(zip(EVK, (Dd.std(0, ddof=1) / np.sqrt(len(Dd))).tolist()))
out["total_seconds"] = time.time() - t_all
json.dump(out, open(f"{_ROOT}/out/cond_{MODE}_N{N}_R{REPS}_s{SEED}.json", "w"), indent=1)
k = "Oth>=5 (actual or more)"
print(f"done {MODE} N={N} R={REPS}: Grn {out['mean_0.0']['Grn>=4 (actual or more)']:.6f}±{out['se_0.0']['Grn>=4 (actual or more)']:.6f} "
      f"Oth {out['mean_0.0'][k]:.6f}±{out['se_0.0'][k]:.6f}  diff(Lab326, lam .5-0) {out['diff_lam05_minus_0']['Lab>=326 (majority)']:.5f}±{out['diff_se_crn']['Lab>=326 (majority)']:.5f}", flush=True)
