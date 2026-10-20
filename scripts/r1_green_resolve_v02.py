"""R1 Green gap, step after the fixed-outer audit (out/r1_fixed_outer_v02.json: regional rule (3,3,3,9,7,3) vs regional
RQMC at 16 fixed outer points, mean diff +5.7e-5 +- 5.6e-5).
  1. Outer check: hybrid (3,3,3,3,3,3) with NEW outer points (2^12 x 8, seed 1), compared with the stored seed-0 run
     (2^11 x 4). The paired rule increment to (3,3,3,9,5,3) (+7.154e-4, paired SE 7e-6, seed 0) is added unchanged.
  2. Reference check: a second independent conditional RQMC reference (N = 8192, R = 32, new seed 20261103), same
     script as the stored one (scripts/cond_mc_v02.py, seed 20261022).
Writes out/r1_green_resolve_v02.json."""
import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import sys, json, time, subprocess, numpy as np
sys.path.insert(0, _ROOT)
OUT = f"{_ROOT}/out/r1_green_resolve_v02.json"
res = json.load(open(OUT)) if _os.path.exists(OUT) else {}
step = sys.argv[1] if len(sys.argv) > 1 else "all"
if step in ("all", "outer"):
    src = open(f"{_ROOT}/scripts/v02_experiments.py").read()
    argv0 = sys.argv
    sys.argv = [sys.argv[0], "4"]
    exec(src[:src.index('out = {"kappa"')])
    sys.argv = argv0
    mp2 = mk(C.POLLS_2024)
    EVK = list(events(QEngine2(mp2, nat_on=False, mu_fixed=np.zeros(4)).run(reps=1, m_log2=0)["dist"]).keys())
    t0 = time.time(); r_ = QEngine2(mp2, reg_ns=(3, 3, 3, 3, 3, 3)).run(m_log2=12, reps=8, seed=1); sec = time.time() - t0
    E = np.array([[events(x)[kk] for kk in EVK] for x in r_["dist_reps"]])
    res["outer_seed1_333333"] = dict(sec=sec, m_log2=12, reps=8, seed=1, mean=dict(zip(EVK, E.mean(0).tolist())),
                                     se=dict(zip(EVK, (E.std(0, ddof=1) / np.sqrt(len(E))).tolist())))
    k = "Grn>=4 (actual or more)"
    print(f"outer seed 1: {sec:.0f}s  Grn {E[:, EVK.index(k)].mean():.6f} +- {E[:, EVK.index(k)].std(ddof=1)/np.sqrt(len(E)):.1e}",
          flush=True)
    json.dump(res, open(OUT, "w"), indent=1)
if step in ("all", "ref"):
    t0 = time.time()
    subprocess.run([sys.executable, f"{_ROOT}/scripts/cond_mc_v02.py", "rqmc", "8192", "32", "20261103"], check=True)
    o = json.load(open(f"{_ROOT}/out/cond_rqmc_N8192_R32_s20261103.json"))
    res["reference_seed20261103"] = dict(sec=time.time() - t0, mean=o["mean_0.0"], se=o["se_0.0"])
    json.dump(res, open(OUT, "w"), indent=1)
print("done", flush=True)
