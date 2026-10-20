"""R1 regional-rule bias with a valid RQMC standard error (reply to the R1 follow-up review of 2026-10-20, item 1).

The earlier audit (r1_fixed_outer_all_v02.py) used 32 points of ONE scrambled Sobol net and took sd/sqrt(32) as the SE,
which is not valid (points of one net are dependent). Here the outer points come from R independent scrambles
(seed [20261102, r]), n points each; the per-replicate mean difference is one observation, and the SE is the sd of the
R replicate means / sqrt(R) (Student-t, R - 1 df).

At each outer point (z, mu): hybrid core with the AUDITED regional rule (3,3,3,9,5,3) (deterministic given the point)
minus regional-only conditional RQMC over the 57 regional factors (2^NREG points x REPS scrambles, seed
[20261102, r, p, k]); the inner RQMC is unbiased, so its noise is part of the replicate variance.
Writes out/r1_fixed_outer_reps_v02_seedfix.json (raw per-point values and the replicate summary). The stored run used by the
tables, out/r1_fixed_outer_reps_v02.json, was written by the pre-fix version (old seeds) and is kept unchanged.

SEED FIX (2026-10-20, R1 follow-up review): the stored run used outer default_rng([SEED, r]) and inner
default_rng([SEED, r, p, k]). numpy's SeedSequence pads short entropy with zeros, so [SEED, r] and [SEED, r, 0, 0] give
the SAME stream: the first outer point's inner draws were not independent of the outer draw. The code now separates the
streams (outer [SEED, 0, r], inner [SEED, 1, r, p, k]); the stored JSON (old seeds) is re-aggregated without each
replicate's first outer point (points[1:], 7-point means) in scripts/r1_note_tables_v02.py.""" 
import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import sys, json, time, numpy as np
from scipy.stats import qmc, norm, t as student_t
sys.path.insert(0, _ROOT)
src = open(f"{_ROOT}/scripts/cond_mc_v02.py").read()
sys.argv = [sys.argv[0], "rqmc", "1024", "1", "0"]
exec(src[:src.index("EVK = None")].replace("cache=True", "cache=False"))
from ukseat.core import exact_chunk_L
R = int(_os.environ.get("R", "16")); NPT = int(_os.environ.get("NPT", "8"))
NREG = int(_os.environ.get("NREG", "12")); REPS = int(_os.environ.get("REPS", "4"))
RULE = (3, 3, 3, 9, 5, 3); SEED = 20261102
OUT = f"{_ROOT}/out/r1_fixed_outer_reps_v02_seedfix.json"  # the stored old-seed run is out/r1_fixed_outer_reps_v02.json; never overwritten
eng = QEngine2(mp, reg_ns=RULE)
rs, rw, rc = eng._reg_arrays()
res = dict(rule=list(RULE), R=R, n_points_per_replicate=NPT, inner=dict(n=2 ** NREG, scrambles=REPS), seed=SEED,
           replicates=[])
EVK = None
t_all = time.time()
for r in range(R):
    Zo = norm.ppf(qmc.Sobol(dz + dm, scramble=True, seed=np.random.default_rng([SEED, 0, r])).random(NPT))
    pts = []
    for p in range(NPT):
        z = Zo[p, :dz]; mu = mp.mu_mean + mp.mu_sd * Zo[p, dz:]
        Dn = (mp.scenario_shift(0.0)[None, :] + (z @ mp.A.T)[None, :])[:, :len(CATS)].copy()
        dd = exact_chunk_L(Dn, eng.base_s, eng.pres_s, eng.L_s, mu[None, :].copy(), eng.roff, rs, rw, rc,
                           mp.pi, mp.tau1, mp.tau2, eng.eps, eng.N)[0]
        EVh = events(dd[0]); EVK = EVK or list(EVh)
        acc_all = []
        eta0 = base0 + (A @ z)[None, :] + np.einsum("ikj,j->ik", L, mu) + shift[0.0][None, :]
        for k in range(REPS):
            U = qmc.Sobol(dr, scramble=True, seed=np.random.default_rng([SEED, 1, r, p, k])).random_base2(NREG)
            W = norm.ppf(np.clip(U, 1e-12, 1 - 1e-12))
            acc = np.zeros((len(CATS), n + 1))
            for row in W:
                eta = eta0.copy(); off = 0
                for idx, Cj in regs:
                    w = row[off: off + Cj.shape[1]]; off += Cj.shape[1]
                    eta[idx] += (Cj @ w)[None, :]
                acc += pb_all(probs(eta))
            acc_all.append([events(acc / len(W))[kk] for kk in EVK])
        mc = np.array(acc_all).mean(0)
        pts.append(dict(hybrid=[EVh[kk] for kk in EVK], mc=mc.tolist()))
    D = np.array([np.array(x["hybrid"]) - np.array(x["mc"]) for x in pts]).mean(0)
    res["replicates"].append(dict(r=r, mean_diff=D.tolist(), points=pts))
    gi = EVK.index("Grn>=4 (actual or more)")
    print(f"replicate {r}: Grn mean diff {D[gi]:+.2e}  ({time.time() - t_all:.0f}s)", flush=True)
    res["events"] = EVK
    json.dump(res, open(OUT, "w"), indent=1)
M = np.array([x["mean_diff"] for x in res["replicates"]])
m, se = M.mean(0), M.std(0, ddof=1) / np.sqrt(len(M)); tq = student_t.ppf(0.975, len(M) - 1)
res["summary"] = {kk: dict(mean=float(m[i]), se=float(se[i]), ci95=[float(m[i] - tq * se[i]), float(m[i] + tq * se[i])])
                  for i, kk in enumerate(EVK)}
res["seconds"] = time.time() - t_all
json.dump(res, open(OUT, "w"), indent=1)
for kk in EVK:
    s = res["summary"][kk]
    print(f"  {kk:30s} {s['mean']:+.2e} +- {s['se']:.1e}  95% [{s['ci95'][0]:+.1e}, {s['ci95'][1]:+.1e}]", flush=True)
print("done", flush=True)
