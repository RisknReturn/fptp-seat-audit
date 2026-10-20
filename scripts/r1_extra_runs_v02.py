"""R1 follow-up: one more independent hybrid outer run (rule (3,3,3,3,3,3), 2^12 x 8, seed 2; per-scramble events saved)
and one more conditional RQMC reference (N = 8192, R = 32, seed 20261104). Motivation: after the valid-SE regional audit,
Lab <= 411 showed a difference whose 95% interval excluded zero; more independent runs show whether it persists.
Pooling rules (set before these runs): hybrid level = equal-weight mean of the 2^12 x 8 runs (seed 1, seed 2);
reference = all scramble replicates of all references, equal weight. Writes out/r1_extra_runs_v02.json."""
import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import sys, json, time, subprocess, numpy as np
sys.path.insert(0, _ROOT)
OUT = f"{_ROOT}/out/r1_extra_runs_v02.json"
res = json.load(open(OUT)) if _os.path.exists(OUT) else {}
step = sys.argv[1] if len(sys.argv) > 1 else "all"
if step in ("all", "outer"):
    src = open(f"{_ROOT}/scripts/v02_experiments.py").read()
    argv0 = sys.argv; sys.argv = [sys.argv[0], "4"]
    exec(src[:src.index('out = {"kappa"')])
    sys.argv = argv0
    mp2 = mk(C.POLLS_2024)
    EVK = list(events(QEngine2(mp2, nat_on=False, mu_fixed=np.zeros(4)).run(reps=1, m_log2=0)["dist"]).keys())
    t0 = time.time(); r_ = QEngine2(mp2, reg_ns=(3, 3, 3, 3, 3, 3)).run(m_log2=12, reps=8, seed=2); sec = time.time() - t0
    E = np.array([[events(x)[kk] for kk in EVK] for x in r_["dist_reps"]])
    res["outer_seed2_333333"] = dict(sec=sec, m_log2=12, reps=8, seed=2, events=EVK, reps_events=E.tolist(),
                                     mean=dict(zip(EVK, E.mean(0).tolist())),
                                     se=dict(zip(EVK, (E.std(0, ddof=1) / np.sqrt(len(E))).tolist())))
    print(f"outer seed 2: {sec:.0f}s", flush=True)
    json.dump(res, open(OUT, "w"), indent=1)
if step in ("all", "ref"):
    t0 = time.time()
    subprocess.run([sys.executable, f"{_ROOT}/scripts/cond_mc_v02.py", "rqmc", "8192", "32", "20261104"], check=True)
    res["reference_seed20261104"] = dict(sec=time.time() - t0, file="out/cond_rqmc_N8192_R32_s20261104.json")
    json.dump(res, open(OUT, "w"), indent=1)
print("done", flush=True)
