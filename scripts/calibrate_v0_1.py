"""Estimate regional-deviation covariance and seat-level noise scale from 2010-2019
(same 2010 boundaries). Uses only pre-2024 information."""
import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import sys, json, numpy as np, pandas as pd
sys.path.insert(0, _ROOT)
from ukseat.swing import shares, solve_delta

h = pd.read_csv(f"{_ROOT}/data/old_bound_hist_v0_1.csv", index_col=0)
YRS = [2010, 2015, 2017, 2019]
# contender categories; 'RefL' = UKIP (2010-17) or Brexit Party (2019)
CATS = ["Con", "Lab", "LD", "RefL", "Grn", "SNP", "PC", "Oth", "R"]

def table(yr):
    t = pd.DataFrame(index=h.index)
    for k in ["Con", "Lab", "LD", "Grn", "SNP", "PC"]:
        t[k] = h[f"{yr}_{k}"]
    t["RefL"] = h[f"{yr}_UKIP"] + h[f"{yr}_Ref"]
    t["Oth"] = h[f"{yr}_OthMax"]
    t["R"] = h[f"{yr}_Oth"] - h[f"{yr}_OthMax"] + h[f"{yr}_Spk"]
    return t[CATS]

V = {y: table(y) for y in YRS}
eng = h.country == "England"
regions = sorted(h.loc[eng, "region"].unique())

# ---------- 1) regional deviations (English regions vs England), log-share space
rows = []
for y0, y1 in zip(YRS[:-1], YRS[1:]):
    E0 = V[y0][eng].sum(); E1 = V[y1][eng].sum()
    sE = np.log((E1 / E1.sum()) / (E0 / E0.sum()))
    for r in regions:
        m = eng & (h.region == r)
        R0 = V[y0][m].sum(); R1 = V[y1][m].sum()
        sR = np.log((R1 / R1.sum()) / (R0 / R0.sum()))
        d = (sR - sE)
        rows.append(dict(t=f"{y0}-{y1}", region=r, **{k: d[k] for k in ["Con", "Lab", "LD", "RefL", "Grn"]}))
dev = pd.DataFrame(rows)
print("Regional log-share deviation sd by transition:")
print(dev.groupby("t")[["Con", "Lab", "LD", "RefL", "Grn"]].std().round(3))
main = dev[["Con", "Lab", "LD"]]
cov3 = np.cov(main.values.T, ddof=1)
print("pooled cov (Con,Lab,LD):\n", cov3.round(5))
print("pooled corr:\n", np.corrcoef(main.values.T).round(2))
# Ref-like and Green: only 2010-15 is free of large standing-pattern artefacts
sd_ref = dev.loc[dev.t == "2010-2015", "RefL"].std()
sd_grn = dev.loc[dev.t == "2010-2015", "Grn"].std()
print("2010-15 regional sd RefL %.3f Grn %.3f" % (sd_ref, sd_grn))

# ---------- 2) seat-level residual scale tau (logistic on top-2 log-ratio residuals)
res_all = []
for y0, y1 in zip(YRS[:-1], YRS[1:]):
    for r in sorted(h.region.unique()):
        m = (h.region == r).values
        v0 = V[y0][m].values; v1 = V[y1][m].values
        s0 = v0 / v0.sum(1, keepdims=True); s1 = v1 / v1.sum(1, keepdims=True)
        with np.errstate(divide="ignore"):
            logb = np.log(s0)
        # regional log share growth: log-ratio residuals only need delta_k - delta_j
        S0 = v0.sum(0) / v0.sum(); S1 = v1.sum(0) / v1.sum()
        with np.errstate(divide="ignore", invalid="ignore"):
            g = np.log(S1 / S0)
        for i in range(m.sum()):
            a = s1[i]
            order = np.argsort(-a)
            k, j = order[0], order[1]
            if CATS[k] in ("R",) or CATS[j] in ("R",):
                continue
            # require both stood at t0 and t1 with >=5% at both
            if min(s0[i, k], s0[i, j], a[k], a[j]) < 0.05:
                continue
            rr = np.log(a[k] / a[j]) - np.log(s0[i, k] / s0[i, j]) - (g[k] - g[j])
            res_all.append(dict(t=f"{y0}-{y1}", region=r, pair="-".join(sorted([CATS[k], CATS[j]])), res=rr))
res = pd.DataFrame(res_all)
def tau_iqr(x):
    q = np.quantile(x, [0.25, 0.75]); return (q[1] - q[0]) / (2 * np.log(3))
print("\nSeat residual (top-2 log-ratio) logistic scale tau by transition (IQR-based):")
for t, g in res.groupby("t"):
    x = g.res.values
    print(f"  {t}: n={len(x)} tau={tau_iqr(x):.3f}  sd={x.std():.3f}  (logistic sd would be {tau_iqr(x)*np.pi/np.sqrt(3):.3f})")
x = res.res.values
tau = tau_iqr(x)
print(f"  pooled: n={len(x)} tau={tau:.3f}")
# tail check: share of |res| > c vs logistic prediction
for c in [1, 2, 3]:
    thr = c * tau * np.pi / np.sqrt(3)
    emp = np.mean(np.abs(x) > thr)
    th = 2 / (1 + np.exp(thr / tau))
    print(f"  P(|res|>{c} sd): empirical {emp:.4f}  logistic {th:.4f}  normal {2*(1-__import__('scipy').stats.norm.cdf(c)):.4f}")

# ---------- 3) two-scale logistic mixture for seat noise (MLE), used by the model
from scipy.optimize import minimize
def logpdf(x, s):
    z = np.abs(x) / s
    return -z - 2 * np.log1p(np.exp(-z)) - np.log(s)
def nll(p):
    pi = 1 / (1 + np.exp(-p[0])); s1 = np.exp(p[1]); s2 = np.exp(p[2])
    return -np.logaddexp(np.log1p(-pi) + logpdf(x, s1), np.log(pi) + logpdf(x, s2)).sum()
best = min((minimize(nll, [np.log(q / (1 - q)), np.log(0.1), np.log(s2)], method="Nelder-Mead",
                     options=dict(maxiter=4000, xatol=1e-8, fatol=1e-8)) for q in (0.05, 0.15, 0.3) for s2 in (0.25, 0.4)),
           key=lambda o: o.fun)
mix_pi = 1 / (1 + np.exp(-best.x[0])); mix_t1, mix_t2 = np.exp(best.x[1:])
one = minimize(lambda p: -logpdf(x, np.exp(p[0])).sum(), [np.log(0.13)])
print("mixture: pi=%.3f tau1=%.3f tau2=%.3f ; single logistic MLE tau=%.3f" % (mix_pi, mix_t1, mix_t2, np.exp(one.x[0])))
from scipy.stats import skew
print("residual skewness %.2f" % skew(x))

out = dict(cov_reg_ConLabLD=cov3.tolist(), sd_reg_RefL=float(sd_ref), sd_reg_Grn=float(sd_grn),
           mix_pi=float(mix_pi), mix_tau1=float(mix_t1), mix_tau2=float(mix_t2), tau_single_mle=float(np.exp(one.x[0])),
           tau_pooled=float(tau),
           tau_by_t={t: float(tau_iqr(g.res.values)) for t, g in res.groupby("t")})
json.dump(out, open(f"{_ROOT}/data/calib_v0_1.json", "w"), indent=1)
dev.to_csv(f"{_ROOT}/data/regional_dev_v0_1.csv", index=False)
res.to_csv(f"{_ROOT}/data/seat_resid_v0_1.csv", index=False)
print(json.dumps(out, indent=1))
