"""Historical realisations (2010->15, 15->17, 17->19, 2010 boundaries) of the seat-type factors used in v0.2.
Residual r_ik = log(v1_ik/v0_ik) - regional log swing of party k. Factor realisations:
  mu_T (anti-Conservative tactical consolidation): Con-held seats, challenger (Lab or LD, larger at t0) vs the other:
        mean(r_chal - r_other)/2      [loading +1 challenger, -1 other]
  mu_C (beyond-proportional Conservative change where Con was strong): OLS slope of r_Con on (v0_Con - regional mean)
  mu_R (Leave -> UKIP/Brexit Party): OLS slope of r_RefL on (Leave-52)/10, seats where it stood both times
  mu_L (Muslim share -> Labour): OLS slope of r_Lab on (Muslim-mean)/10
Census/Leave are on 2024 boundaries; mapped back to 2010 boundaries with population overlaps (approximation)."""
import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import sqlite3, json, numpy as np, pandas as pd
h = pd.read_csv(f"{_ROOT}/data/old_bound_hist_v0_1.csv", index_col=0)
cov = pd.read_csv(f"{_ROOT}/data/gb2024_v02_covariates.csv", index_col=0)
c = sqlite3.connect(_os.environ.get("UKSEAT_DB", _os.path.join(_ROOT, "..", "data_raw", "psephology.db")))
ov = pd.read_sql("select from_constituency_area_id f, to_constituency_area_id t, from_constituency_population fp from constituency_area_overlaps", c)
ov = ov[ov.f.isin(h.index) & ov.t.isin(cov.index)]
mapped = {}
for v in ["leave", "muslim", "degree"]:
    x = ov.merge(cov[[v]], left_on="t", right_index=True)
    mapped[v] = (x.fp * x[v]).groupby(x.f).sum() / x.fp.groupby(x.f).sum()
h = h.join(pd.DataFrame(mapped))
YRS = [2010, 2015, 2017, 2019]
def shares(y):
    t = pd.DataFrame(index=h.index)
    for k in ["Con", "Lab", "LD", "Grn"]:
        t[k] = h[f"{y}_{k}"]
    t["RefL"] = h[f"{y}_UKIP"] + h[f"{y}_Ref"]
    t["tot"] = h[[c_ for c_ in h.columns if c_.startswith(f"{y}_") and not c_.endswith("OthMax")]].sum(axis=1)
    return t
res = {}
for y0, y1 in zip(YRS[:-1], YRS[1:]):
    a, b = shares(y0), shares(y1)
    v0 = a[["Con", "Lab", "LD", "RefL", "Grn"]].div(a.tot, axis=0); v1 = b[["Con", "Lab", "LD", "RefL", "Grn"]].div(b.tot, axis=0)
    r = pd.DataFrame(index=h.index)
    for k in v0:
        ok = (v0[k] > 0.02) & (v1[k] > 0.02)
        R0 = a[k].groupby(h.region).sum() / a.tot.groupby(h.region).sum(); R1 = b[k].groupby(h.region).sum() / b.tot.groupby(h.region).sum()
        g = np.log(R1 / R0)
        r[k] = np.where(ok, np.log(v1[k] / v0[k]) - h.region.map(g), np.nan)
    eng_wal = h.country != "Scotland"
    # mu_T
    conheld = (v0.idxmax(axis=1) == "Con") & eng_wal
    chal_is_lab = v0.Lab >= v0.LD
    d = np.where(chal_is_lab, r.Lab - r.LD, r.LD - r.Lab)
    muT = np.nanmean(d[conheld]) / 2
    # mu_C
    x = (v0.Con - v0.Con.groupby(h.region).transform("mean")); okc = r.Con.notna() & eng_wal
    muC = np.polyfit(x[okc], r.Con[okc], 1)[0]
    # mu_R
    xr = (h.leave - 52) / 10; okr = r.RefL.notna() & xr.notna()
    muR = np.polyfit(xr[okr], r.RefL[okr], 1)[0] if okr.sum() > 50 else np.nan
    # mu_L (Muslim -> Labour)
    xm = (h.muslim - h.muslim.mean()) / 10; okm = r.Lab.notna() & xm.notna()
    muL = np.polyfit(xm[okm], r.Lab[okm], 1)[0]
    res[f"{y0}-{y1}"] = dict(mu_T=float(muT), mu_C=float(muC), mu_R=float(muR), mu_L=float(muL),
                            n_conheld=int(conheld.sum()), n_ref=int(okr.sum()))
    print(f"{y0}-{y1}", {k: round(v, 3) if isinstance(v, float) else v for k, v in res[f'{y0}-{y1}'].items()})
summ = {}
for k in ["mu_T", "mu_C", "mu_R", "mu_L"]:
    v = np.array([res[t][k] for t in res]); v = v[np.isfinite(v)]
    summ[k] = dict(mean=float(v.mean()), sd=float(v.std(ddof=1)) if len(v) > 1 else float("nan"), values=v.tolist())
    print(k, "mean %.3f sd %.3f" % (summ[k]["mean"], summ[k]["sd"]))
json.dump(dict(by_transition=res, summary=summ), open(f"{_ROOT}/data/calib_mu_v02.json", "w"), indent=1)
