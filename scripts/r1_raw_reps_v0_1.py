"""Raw replicate values for the two N7 inputs whose stored files kept only means and SEs (review of 2026-10-21, item 2).
Both are re-runs with the original seeds; the engine is deterministic given the seed, so they must reproduce the stored
values exactly.
  level:      hybrid rule (3,3,3,3,3,3), 2^12 x 8, seed 1 (r1_green_resolve_v02.py, step "outer"); checked against the
              stored mean and SE in r1_green_resolve_v02.json.
  increment:  rule (3,3,3,9,5,3), 2^11 x 4, seed 0, lambda = 0, paired with the stored (3,3,3,3,3,3) replicates at the
              same points (reg_convergence_v02b.json); checked against the stored increment and its paired SE in
              r1_followup_v02.json.
usage: python r1_raw_reps_v0_1.py [level|increment|all]   (default all; a partial run keeps the other entry)
Writes out/r1_raw_reps_v0_1.json. Does not modify any other stored file."""
import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import sys, json, time, numpy as np
sys.path.insert(0, _ROOT)
step = sys.argv[1] if len(sys.argv) > 1 else "all"
src = open(f"{_ROOT}/scripts/v02_experiments.py").read()
sys.argv = [sys.argv[0], "4"]
exec(src[:src.index('out = {"kappa"')])
mp2 = mk(C.POLLS_2024)
EVK = list(events(QEngine2(mp2, nat_on=False, mu_fixed=np.zeros(4)).run(reps=1, m_log2=0)["dist"]).keys())
O = f"{_ROOT}/out"
OUT = f"{O}/r1_raw_reps_v0_1.json"
res = json.load(open(OUT)) if _os.path.exists(OUT) else {}
res["events"] = EVK


def reps_of(rule, **kw):
    t0 = time.time(); r_ = QEngine2(mp2, reg_ns=rule).run(**kw)
    return np.array([[events(x)[k] for k in EVK] for x in r_["dist_reps"]]), time.time() - t0


if step in ("all", "level"):
    st = json.load(open(f"{O}/r1_green_resolve_v02.json"))["outer_seed1_333333"]
    E, sec = reps_of((3, 3, 3, 3, 3, 3), m_log2=12, reps=8, seed=1)
    dmean = np.abs(E.mean(0) - np.array([st["mean"][k] for k in EVK])).max()
    dse = np.abs(E.std(0, ddof=1) / np.sqrt(len(E)) - np.array([st["se"][k] for k in EVK])).max()
    res["level_seed1_333333"] = dict(rule=[3, 3, 3, 3, 3, 3], m_log2=12, reps=8, seed=1, sec=sec, reps_events=E.tolist(),
                                     max_abs_diff_mean_vs_stored=float(dmean), max_abs_diff_se_vs_stored=float(dse))
    print(f"level, seed 1: {sec:.0f}s, max |mean - stored| {dmean:.1e}, max |SE - stored| {dse:.1e}", flush=True)
    json.dump(res, open(OUT, "w"), indent=1)

if step in ("all", "increment"):
    rg = json.load(open(f"{O}/reg_convergence_v02b.json"))
    E333 = np.array(rg["paired"]["(3, 3, 3, 3, 3, 3)"]["reps"])
    Eaud, sec = reps_of((3, 3, 3, 9, 5, 3), lam=0.0, m_log2=11, reps=4, seed=0)
    D = Eaud - E333
    fu = json.load(open(f"{O}/r1_followup_v02.json"))["green_axes"]["(3, 3, 3, 9, 5, 3)"]
    dm = np.abs(D.mean(0) - np.array([fu["diff_vs_333333"][k] for k in EVK])).max()
    ds = np.abs(D.std(0, ddof=1) / np.sqrt(len(D)) - np.array([fu["diff_se"][k] for k in EVK])).max()
    res["increment_seed0_333333_to_333953"] = dict(m_log2=11, reps=4, seed=0, sec_333953=sec, reps_333333=E333.tolist(),
                                                   reps_333953=Eaud.tolist(), reps_diff=D.tolist(),
                                                   max_abs_diff_mean_vs_stored=float(dm), max_abs_diff_se_vs_stored=float(ds))
    print(f"increment, seed 0: {sec:.0f}s, max |mean - stored| {dm:.1e}, max |SE - stored| {ds:.1e}", flush=True)
    json.dump(res, open(OUT, "w"), indent=1)
print("done", flush=True)
