"""R2 review 1.1: the 'realised mu' numbers are descriptive statistics h(mu), not generative coefficients.

1. Synthetic recovery: apply the descriptive estimator (same logic as realised_mu_2024.py) to noiseless votes generated
   by the model itself (same candidates, baseline, loadings, exact re-raking) at mu = 0 and at mu = the statistics.
2. Jacobian of h at 0.
3. Inversion: solve h(mu*) = h_obs (statistics of the actual 2024 votes, computed on the same seat set) by Newton
   steps with a central-difference Jacobian. mu* is the coefficient vector whose generated votes reproduce the observed
   statistics (method-of-moments), still a 2024 development diagnostic, not a forecast.
4. Oracle seats (actual GB/Scotland/Wales vote shares, national error off, mu fixed and re-raked) at mu = 0,
   mu = h_obs (the v0.2 practice) and mu = mu*, with regional rule (3,3,3,3,3,3); Brier and log score of seat winners.
Writes out/mu_inversion_v02.json. KAPPA = 4 as in v0.2.
"""
import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import sys, json, time, numpy as np, pandas as pd
sys.path.insert(0, _ROOT)
src = open(f"{_ROOT}/scripts/v02_experiments.py").read()
sys.argv = [sys.argv[0], "4"]
exec(src[:src.index('out = {"kappa"')])
from ukseat.model import ALL
from ukseat.swing import shares
from ukseat.quad import tensor_gh

m = mk(POLLS_ACT)
s = m.df
cov = m.cov
BASE_DBAR = m.dbar.copy()
NAMES = ["Con", "Lab", "LD", "RefL", "Grn"]


def tab19():
    a = pd.DataFrame(index=s.index)
    for k in ["Con", "Lab", "LD", "Grn"]:
        a[k] = s[f"n19_{k}"]
    a["RefL"] = s.n19_Ref + s.n19_UKIP
    den = s[[c for c in s.columns if c.startswith("n19_") and not c.endswith("OthMax")]].sum(axis=1)
    return a, den


A19, DEN19 = tab19()
V0 = A19.div(DEN19, axis=0)


def stats_from_shares(v1):
    """Descriptive statistics of realised_mu_2024.py; v1 = 2024 shares (DataFrame with NAMES columns) on the model seats."""
    r = pd.DataFrame(index=s.index)
    for k in NAMES:
        ok = (V0[k] > 0.02) & (v1[k] > 0.02)
        R0 = A19[k].groupby(s.region).sum() / DEN19.groupby(s.region).sum()
        R1 = (v1[k] * DEN19).groupby(s.region).sum() / DEN19.groupby(s.region).sum()
        with np.errstate(divide="ignore", invalid="ignore"):
            r[k] = np.where(ok, np.log(v1[k] / V0[k]) - s.region.map(np.log(R1 / R0)), np.nan)
    ew = s.country != "Scotland"
    held = (V0.idxmax(axis=1) == "Con") & ew
    d = np.where(V0.Lab >= V0.LD, r.Lab - r.LD, r.LD - r.Lab)
    out = [np.nanmean(d[held]) / 2]
    for col, x, mask in (("Con", V0.Con - V0.Con.groupby(s.region).transform("mean"), ew),
                         ("RefL", (cov.leave - 52) / 10, pd.Series(True, index=s.index)),
                         ("Lab", (cov.muslim - cov.muslim.mean()) / 10, pd.Series(True, index=s.index))):
        ok = r[col].notna() & x.notna() & mask
        out.append(np.polyfit(x[ok], r.loc[ok, col], 1)[0])
    return np.array(out)


def gen_shares(mu):
    m.dbar = BASE_DBAR.copy()
    if np.any(mu != 0):
        m.dbar = m.rerake(mu)
    V = shares(m.logb, m.mean_delta() + m.L @ mu)
    return pd.DataFrame({k: V[:, IDX["Ref" if k == "RefL" else k]] for k in NAMES}, index=s.index)


def h(mu):
    return stats_from_shares(gen_shares(np.asarray(mu, float)))


tot24 = s[[c for c in s.columns if c.startswith("a24_") and not c.endswith("OthMax")]].sum(axis=1)
v24 = pd.DataFrame({k: (s[f"a24_{k}"] if k != "RefL" else s.a24_Ref + s.a24_UKIP) / tot24 for k in NAMES}, index=s.index)
h_obs = stats_from_shares(v24)
res = dict(factors=FACTORS, h_obs=h_obs.tolist(), realised_json=MU_REAL.tolist(), kappa=KAPPA)
res["synthetic_recovery"] = {"mu=0": h(np.zeros(4)).tolist(), "mu=h_obs": h(h_obs).tolist()}
print("h_obs", h_obs.round(4), "| h(0)", np.round(res["synthetic_recovery"]["mu=0"], 4),
      "| h(h_obs)", np.round(res["synthetic_recovery"]["mu=h_obs"], 4), flush=True)


def jac(mu, step=1e-4):
    J = np.zeros((4, 4))
    for j in range(4):
        e = np.zeros(4); e[j] = step
        J[:, j] = (h(mu + e) - h(mu - e)) / (2 * step)
    return J


J0 = jac(np.zeros(4))
res["jacobian_at_0"] = J0.tolist()
mu = h_obs.copy()
hist = []
for it in range(12):
    g = h(mu) - h_obs
    hist.append(dict(it=it, mu=mu.tolist(), resid=g.tolist(), max_abs=float(np.abs(g).max())))
    print(f"newton {it}: mu {mu.round(4)} |h-h_obs| {np.abs(g).max():.2e}", flush=True)
    if np.abs(g).max() < 1e-9:
        break
    J = jac(mu)
    mu = mu - np.linalg.solve(J, g)
Jstar = jac(mu)
res["inversion"] = dict(mu_star=mu.tolist(), history=hist, jacobian_at_mu_star=Jstar.tolist(),
                        cond_jacobian=float(np.linalg.cond(Jstar)))
m.dbar = BASE_DBAR.copy()

# oracle seats and seat-level scores under regional rule (3,3,3,3,3,3)
K = len(CATS)
pres = np.isfinite(m.logb[:, :K])
regs = []
for j, Cj in enumerate(m.C):
    ix = np.where(m.reg == j)[0]
    z, w = tensor_gh((3, 3, 3, 3, 3, 3)[:Cj.shape[1]])
    regs.append((ix, z @ Cj[:K].T, w))
act = s[[f"a24_{k}" for k in ["Con", "Lab", "LD", "Ref", "Grn", "SNP", "PC", "OthMax"]]].values.argmax(1)


def oracle(mu):
    m.dbar = BASE_DBAR.copy()
    if np.any(mu != 0):
        m.dbar = m.rerake(mu)
    eta = m.logb[:, :K] + m.mean_delta()[:, :K] + m.L[:, :K] @ mu
    P = np.zeros((m.n, K))
    for ix, rs, rw in regs:
        P[ix] = np.tensordot(rw, m.win_probs(eta[ix][None, :, :] + rs[:, None, :], pres[ix][None, :, :]), axes=1)
    Y = np.zeros_like(P); Y[np.arange(m.n), act] = 1
    gb = shares(m.logb, m.mean_delta() + m.L @ mu)
    out = dict(mean_seats=dict(zip(CATS, P.sum(0).round(3).tolist())), brier=float(((P - Y) ** 2).sum(1).mean()),
               logscore=float(-np.log(np.maximum(P[np.arange(m.n), act], 1e-300)).mean()),
               gb_share_pct=dict(zip(ALL, ((m.w @ gb) / m.w.sum() * 100).round(4).tolist())))
    m.dbar = BASE_DBAR.copy()
    return out


res["oracle"] = {}
for tag, mv in (("mu=0", np.zeros(4)), ("mu=h_obs (v0.2 practice)", h_obs), ("mu=mu_star (inverted)", mu)):
    res["oracle"][tag] = oracle(np.asarray(mv, float))
    o = res["oracle"][tag]
    print(f"{tag:28s} Lab {o['mean_seats']['Lab']:.1f} Con {o['mean_seats']['Con']:.1f} LD {o['mean_seats']['LD']:.1f} "
          f"Ref {o['mean_seats']['Ref']:.1f} SNP {o['mean_seats']['SNP']:.1f} | Brier {o['brier']:.4f} log {o['logscore']:.4f}", flush=True)
res["actual_seats"] = ACT
res["note"] = "2024 development diagnostic; mu* is fitted to 2024 outcomes and is not a forecast"
json.dump(res, open(f"{_ROOT}/out/mu_inversion_v02.json", "w"), indent=1)
print("done", flush=True)
