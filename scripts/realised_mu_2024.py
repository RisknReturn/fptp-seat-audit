"""Ex-post (oracle) realisation of the same factors for 2019 notional -> 2024 actual, on 2024 boundaries."""
import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import json, numpy as np, pandas as pd
s = pd.read_csv(f"{_ROOT}/data/gb2024_seats_v0_1.csv", index_col=0)
cov = pd.read_csv(f"{_ROOT}/data/gb2024_v02_covariates.csv", index_col=0)
s = s.join(cov[["leave", "muslim"]])
s = s[s.n19_Spk == 0]
def tab(p):
    t = pd.DataFrame(index=s.index)
    for k in ["Con", "Lab", "LD", "Grn"]:
        t[k] = s[f"{p}_{k}"]
    t["RefL"] = s[f"{p}_Ref"] + s[f"{p}_UKIP"]
    t["tot"] = s[[c for c in s.columns if c.startswith(p + "_") and not c.endswith("OthMax")]].sum(axis=1)
    return t
a, b = tab("n19"), tab("a24")
v0 = a[["Con", "Lab", "LD", "RefL", "Grn"]].div(a.tot, axis=0); v1 = b[["Con", "Lab", "LD", "RefL", "Grn"]].div(b.tot, axis=0)
r = pd.DataFrame(index=s.index)
for k in v0:
    ok = (v0[k] > 0.02) & (v1[k] > 0.02)
    R0 = a[k].groupby(s.region).sum() / a.tot.groupby(s.region).sum(); R1 = b[k].groupby(s.region).sum() / b.tot.groupby(s.region).sum()
    r[k] = np.where(ok, np.log(v1[k] / v0[k]) - s.region.map(np.log(R1 / R0)), np.nan)
ew = s.country != "Scotland"
conheld = (v0.idxmax(axis=1) == "Con") & ew
d = np.where(v0.Lab >= v0.LD, r.Lab - r.LD, r.LD - r.Lab)
muT = np.nanmean(d[conheld]) / 2
x = v0.Con - v0.Con.groupby(s.region).transform("mean"); okc = r.Con.notna() & ew
muC = np.polyfit(x[okc], r.Con[okc], 1)[0]
xr = (s.leave - 52) / 10; okr = r.RefL.notna()
muR = np.polyfit(xr[okr], r.RefL[okr], 1)[0]
xm = (s.muslim - s.muslim.mean()) / 10; okm = r.Lab.notna()
muL = np.polyfit(xm[okm], r.Lab[okm], 1)[0]
hist = json.load(open(f"{_ROOT}/data/calib_mu_v02.json"))["summary"]
out = {}
for k, v in dict(mu_T=muT, mu_C=muC, mu_R=muR, mu_L=muL).items():
    h = hist[k]; z = (v - h["mean"]) / h["sd"]
    out[k] = dict(realised_2024=float(v), hist_mean=h["mean"], hist_sd=h["sd"], z=float(z))
    print(f"{k}: 2024 realised {v:+.3f} | 2010-19 mean {h['mean']:+.3f} sd {h['sd']:.3f} | z = {z:+.1f}")
json.dump(out, open(f"{_ROOT}/data/realised_mu_2024.json", "w"), indent=1)
