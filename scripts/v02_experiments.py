"""v0.2 experiments: oracle with mu (0 / realised / ex-ante), poll-based ex-ante, lambda x mu_T grid, seat calibration."""
import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import sys, json, time, numpy as np, pandas as pd
sys.path.insert(0, _ROOT)
from ukseat.model2 import SeatModel2, FACTORS
from ukseat.setup2024 import reg_cfg
from ukseat import config_2024 as C
from ukseat.engine2 import QEngine2
from ukseat.report import events, summarise_dist
from ukseat.model import CATS, IDX

ROOT = _ROOT
POLLS_ACT = {"GB": C.ACTUAL_GB_2024,
             "Scotland": {"SNP": 30.0, "Lab": 35.3, "Con": 12.7, "LD": 9.7, "Grn": 3.8, "Ref": 7.0},
             "Wales": {"Lab": 37.0, "Con": 18.2, "PC": 14.8, "LD": 6.5, "Grn": 4.7, "Ref": 16.9}}
REAL = json.load(open(f"{ROOT}/data/realised_mu_2024.json"))
MU_REAL = np.array([REAL[f]["realised_2024"] for f in FACTORS])
KAPPA = float(sys.argv[1]) if len(sys.argv) > 1 else 4.0
NEUTRAL = (sys.argv[2] != "raw") if len(sys.argv) > 2 else True
ACT = C.ACTUAL_SEATS_GB_2024

def mk(polls):
    return SeatModel2(f"{ROOT}/data/gb2024_seats_v0_1.csv", f"{ROOT}/data/calib_v0_1.json", polls, C.SIGMA_NAT, C.CORR_NAT,
                      reg_cfg(), e_ref=C.E_REF_2015, cov_csv=f"{ROOT}/data/gb2024_v02_covariates.csv", kappa=KAPPA, neutral=NEUTRAL)

def summ(r, tag):
    d = r["dist"]; s = summarise_dist(d)
    pit = {k: float(d[CATS.index(k), : ACT[k] + 1].sum()) for k in CATS if k in ACT}
    se = {}
    if r["dist_reps"].shape[0] > 1:
        evr = [events(x) for x in r["dist_reps"]]
        se = {k: float(np.std([e[k] for e in evr], ddof=1) / np.sqrt(len(evr))) for k in evr[0]}
    line = " ".join(f"{k} {s[k]['mean']:.0f}({s[k]['q05']}-{s[k]['q95']})" for k in ["Lab", "Con", "LD", "Ref", "SNP", "Oth"])
    print(f"{tag:42s} {line}  PIT Lab {pit['Lab']:.2f} Con {pit['Con']:.2f} LD {pit['LD']:.2f}  ({r['seconds']:.0f}s)", flush=True)
    return dict(summary=s, events=events(d), events_se=se, pit=pit, dist=d.tolist(), seconds=r["seconds"])

def seat_marginals(model, lam, nat_on, mu_mode, mu_mean=None, mu_sd=None, ndraw=3000, seed=3):
    """E[p_ik] over factors by Monte Carlo on factors only (closed-form seat probabilities: no seat-level noise)."""
    rng = np.random.default_rng(seed)
    K = len(CATS); pres = np.isfinite(model.logb[:, :K])
    base = np.where(pres, model.logb[:, :K], 0.0) + model.mean_delta()[:, :K] + model.scenario_shift(lam)[:K]
    mm = model.mu_mean if mu_mean is None else mu_mean; ss = model.mu_sd if mu_sd is None else mu_sd
    P = np.zeros((model.n, K))
    for _ in range(ndraw):
        eta = base.copy()
        if nat_on:
            eta += (model.A @ rng.standard_normal(model.A.shape[1]))[:K]
        for j, C_ in enumerate(model.C):
            idx = model.reg == j
            eta[idx] += (C_ @ rng.standard_normal(C_.shape[1]))[:K]
        mu = mu_mode if isinstance(mu_mode, np.ndarray) else (mm + ss * rng.standard_normal(len(mm)) if mu_mode == "random" else np.zeros(len(mm)))
        eta += np.einsum("ikj,j->ik", model.L[:, :K, :], mu)
        P += model.win_probs(eta, pres)
    return P / ndraw

def calib(model, P):
    act = model.df[[f"a24_{k}" for k in ["Con", "Lab", "LD", "Ref", "Grn", "SNP", "PC", "OthMax"]]].values.argmax(1)
    y = np.zeros_like(P); y[np.arange(len(act)), act] = 1
    return dict(brier=float(((P - y) ** 2).sum(1).mean()), hit=int((P.argmax(1) == act).sum()), exp_hit=float(P.max(1).sum()),
                logscore=float(-np.log(np.maximum(P[np.arange(len(act)), act], 1e-300)).mean()))

class Raked:
    """with Raked(model, mu): exact block re-raking at fixed mu for the duration of a run (peer review 1.2)."""
    def __init__(self, model, mu): self.m, self.mu = model, mu
    def __enter__(self):
        self.old = self.m.dbar.copy()
        if NEUTRAL and self.mu is not None and np.any(np.asarray(self.mu) != 0):
            self.m.dbar = self.m.rerake(self.mu)
    def __exit__(self, *a): self.m.dbar = self.old

def runfix(model, mf, **kw):
    nat = kw.pop("nat_on", True)
    with Raked(model, mf):
        return QEngine2(model, nat_on=nat, mu_fixed=mf).run(**kw)

out = {"kappa": KAPPA, "neutral": NEUTRAL, "mu_realised": MU_REAL.tolist()}
mo = mk(POLLS_ACT)
print("ex-ante mu mean", mo.mu_mean.round(3), "sd", mo.mu_sd.round(3), "| realised", MU_REAL.round(3), flush=True)
E = lambda **k: QEngine2(mo, nat_on=False, **k)
out["oracle_mu0"] = summ(runfix(mo, np.zeros(4), nat_on=False, reps=1, m_log2=0), "oracle, mu=0 (nominations only)")
out["oracle_mu_real"] = summ(runfix(mo, MU_REAL, nat_on=False, reps=1, m_log2=0), "oracle, mu=realised 2024 (ex post)")
out["oracle_mu_exante"] = summ(E().run(m_log2=9, reps=4), "oracle, mu ~ N(2010-19)")
out["oracle_mu_exante_x3"] = summ(E().run(m_log2=9, reps=4, mu_sd=3 * mo.mu_sd), "oracle, mu ~ N(2010-19, 3x sd)")
for j, f in enumerate(FACTORS):  # one realised factor at a time (others at 0) -> contribution
    mf = np.zeros(4); mf[j] = MU_REAL[j]
    out[f"oracle_only_{f}"] = summ(runfix(mo, mf, nat_on=False, reps=1, m_log2=0), f"oracle, only {f} realised")
mp = mk(C.POLLS_2024)
out["poll_mu_exante"] = summ(QEngine2(mp).run(m_log2=11, reps=4), "final polls, lam=0, mu ~ N(2010-19)")
out["poll_mu_real"] = summ(runfix(mp, MU_REAL, m_log2=11, reps=4), "final polls, lam=0, mu=realised")
grid = {}
for lam in (0.0, 0.5, 1.0):
    for muT in (mp.mu_mean[0], 0.10, MU_REAL[0]):
        mf = mp.mu_mean.copy(); mf[0] = muT
        r = runfix(mp, mf, lam=lam, m_log2=10, reps=2)
        s = summarise_dist(r["dist"]); ev = events(r["dist"])
        grid[f"lam={lam}|muT={muT:.3f}"] = dict(Lab=s["Lab"]["mean"], Con=s["Con"]["mean"], LD=s["LD"]["mean"],
                                                 P_Lab_le_411=ev["Lab<=411 (actual or fewer)"], P_Con_ge_121=ev["Con>=121 (actual or more)"])
        print(f"grid lam={lam} muT={muT:+.3f}: Lab {s['Lab']['mean']:.0f} Con {s['Con']['mean']:.0f} LD {s['LD']['mean']:.0f}", flush=True)
out["grid_lam_muT"] = grid
cal = {}
for tag, model, lam, nat, mode in (("oracle_mu0", mo, 0, False, "zero"), ("oracle_mu_real", mo, 0, False, MU_REAL),
                                   ("oracle_mu_exante", mo, 0, False, "random"), ("poll_mu_exante", mp, 0, True, "random"),
                                   ("poll_mu_real", mp, 0, True, MU_REAL)):
    with Raked(model, mode if isinstance(mode, np.ndarray) else None):
        P = seat_marginals(model, lam, nat, mode, ndraw=1500)
    cal[tag] = calib(model, P)
    print("calib", tag, cal[tag], flush=True)
out["calibration"] = cal
json.dump(out, open(f"{ROOT}/out/v02_experiments_k{KAPPA:g}{'' if NEUTRAL else '_raw'}.json", "w"))
print("done")
