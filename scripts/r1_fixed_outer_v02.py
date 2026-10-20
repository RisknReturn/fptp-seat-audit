"""R1 follow-up: is the remaining Green gap (refined hybrid vs conditional RQMC reference, about +4e-4) in the
regional or in the outer integration? At a few FIXED outer points (z, mu), compare
  hybrid regional Gauss-Hermite with the refined rule (3,3,3,9,7,3)   vs
  regional RQMC over the 57 regional factors only (no regional quadrature), same model and PB kernel.
If they agree at fixed outer points, the gap comes from the outer (z, mu) integration of the hybrid or from the reference.
Writes out/r1_fixed_outer_v02.json.
HISTORY ONLY (2026-10-20): not used in any table. This script seeds the outer net with default_rng([20261031, 0]) and the
inner scrambles with default_rng([20261031, p, r]); numpy's SeedSequence zero-pads short entropy, so the inner stream at
p = 0, r = 0 equals the outer stream. Its single-net SE (sd / sqrt(n)) is also not a valid RQMC SE. See
scripts/r1_fixed_outer_reps_v02.py (seed-separated version) and scripts/r1_note_tables_v02.py."""
import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import sys, json, time, numpy as np
from scipy.stats import qmc, norm
sys.path.insert(0, _ROOT)
src = open(f"{_ROOT}/scripts/cond_mc_v02.py").read()
sys.argv = [sys.argv[0], "rqmc", "1024", "1", "0"]
exec(src[:src.index("EVK = None")].replace("cache=True", "cache=False"))
from ukseat.core import exact_chunk_L
NPTS = int(_os.environ.get("NPTS", "6")); NREG = int(_os.environ.get("NREG", "13")); REPS = int(_os.environ.get("REPS", "8"))
RULE = (3, 3, 3, 9, 7, 3)
eng = QEngine2(mp, reg_ns=RULE)
rs, rw, rc = eng._reg_arrays()
gi = CATS.index("Grn")
Zo = norm.ppf(qmc.Sobol(dz + dm, scramble=True, seed=np.random.default_rng([20261031, 0])).random(NPTS))
res = dict(rule=list(RULE), npts=NPTS, nreg=2 ** NREG, reps=REPS, points=[])


def grn_ge4(dist):
    return float(dist[gi, 4:].sum())


for p in range(NPTS):
    z = Zo[p, :dz]; mu = mp.mu_mean + mp.mu_sd * Zo[p, dz:]
    Dn = (mp.scenario_shift(0.0)[None, :] + (z @ mp.A.T)[None, :])[:, :len(CATS)].copy()
    t0 = time.time()
    dd, _, _, _, _ = exact_chunk_L(Dn, eng.base_s, eng.pres_s, eng.L_s, mu[None, :].copy(), eng.roff, rs, rw, rc,
                                   mp.pi, mp.tau1, mp.tau2, eng.eps, eng.N)
    h = grn_ge4(dd[0]); th = time.time() - t0
    t0 = time.time(); vals = []
    for r in range(REPS):
        U = qmc.Sobol(dr, scramble=True, seed=np.random.default_rng([20261031, p, r])).random_base2(NREG)
        W = norm.ppf(np.clip(U, 1e-12, 1 - 1e-12))
        acc = np.zeros((len(CATS), n + 1))
        eta0 = base0 + (A @ z)[None, :] + np.einsum("ikj,j->ik", L, mu) + shift[0.0][None, :]
        for row in W:
            eta = eta0.copy(); off = 0
            for idx, Cj in regs:
                w = row[off: off + Cj.shape[1]]; off += Cj.shape[1]
                eta[idx] += (Cj @ w)[None, :]
            acc += pb_all(probs(eta))
        vals.append(grn_ge4(acc / len(W)))
    vals = np.array(vals); tm = time.time() - t0
    rec = dict(point=p, hybrid=h, regional_rqmc=float(vals.mean()), regional_rqmc_se=float(vals.std(ddof=1) / np.sqrt(REPS)),
               diff=float(h - vals.mean()), sec_hybrid=th, sec_mc=tm)
    res["points"].append(rec)
    print(f"point {p}: hybrid {h:.6f}  regional RQMC {vals.mean():.6f} +- {rec['regional_rqmc_se']:.1e}  diff {rec['diff']:+.2e}"
          f"  ({th:.1f}s / {tm:.0f}s)", flush=True)
    json.dump(res, open(f"{_ROOT}/out/r1_fixed_outer_v02.json", "w"), indent=1)
d = np.array([x["diff"] for x in res["points"]]); se = np.array([x["regional_rqmc_se"] for x in res["points"]])
res["mean_diff"] = float(d.mean()); res["mean_diff_se_mc_only"] = float(np.sqrt((se ** 2).sum()) / len(se))
res["mean_diff_se_points"] = float(d.std(ddof=1) / np.sqrt(len(d))) if len(d) > 1 else None   # includes outer-point spread
res["n_diff_gt_2se"] = int(np.sum(np.abs(d) > 2 * se))
print("mean diff", res["mean_diff"], "+- (MC)", res["mean_diff_se_mc_only"], "+- (points)", res["mean_diff_se_points"],
      "| points with |diff| > 2 SE:", res["n_diff_gt_2se"], flush=True)
json.dump(res, open(f"{_ROOT}/out/r1_fixed_outer_v02.json", "w"), indent=1)
