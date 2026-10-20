"""2024 backtest: final polls -> seat distribution; lambda bias scenarios (2015-size miss multiples);
oracle run with actual national shares; seat-level calibration; variance decomposition."""
import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import sys, json, time, numpy as np
sys.path.insert(0, _ROOT)
from ukseat.setup2024 import make_model
from ukseat.engine import FastEngine, ExactEngine
from ukseat.report import events, summarise_dist
from ukseat.model import CATS
from ukseat import config_2024 as C

REF_NAT = [("trap", 0.25, 7.5), ("gh", 9), ("gh", 5), ("gh", 3), ("gh", 3)]
SCAN_NAT = [("trap", 0.35, 7.5), ("gh", 5), ("gh", 3), ("gh", 5), ("gh", 3)]
REG = sys.argv[1] if len(sys.argv) > 1 else "auto"
REG = REG if REG == "auto" else tuple(int(x) for x in REG.split(","))
OUT = f"{_ROOT}/out"
ACT = C.ACTUAL_SEATS_GB_2024
res = {}

def pack(r):
    sm = summarise_dist(r["dist"])
    sd = np.sqrt(r["var_nat"] + r["var_reg"] + r["var_seat"])
    pit = {k: float(r["dist"][i, : ACT[k] + 1].sum()) for i, k in enumerate(CATS) if k in ACT}
    return dict(summary=sm, events=events(r["dist"]), pit_le_actual=pit, dist=r["dist"].tolist(),
                var_share={k: dict(nat=float(r["var_nat"][i] / sd[i] ** 2), reg=float(r["var_reg"][i] / sd[i] ** 2),
                                   seat=float(r["var_seat"][i] / sd[i] ** 2)) for i, k in enumerate(CATS)},
                seconds=r["seconds"], nodes=r["nodes_outer"], tv_bound=r["tv_bound"])

m = make_model()
# ---- reference run at lambda = 0: taken from scripts/ref_exact2_v0_1.py (all dimensions refined)
r2 = json.load(open(f"{OUT}/ref_exact2_v0_1.json"))
d2 = np.array(r2["dist"]); sd2 = np.sqrt(np.array(r2["var"]["nat"]) + np.array(r2["var"]["reg"]) + np.array(r2["var"]["seat"]))
res["lam0_ref"] = dict(summary=r2["summ"], events=r2["ev"], dist=r2["dist"], seconds=r2["sec"], nodes=r2["nodes"], tv_bound=r2["tv"],
                       pit_le_actual={k: float(d2[i, : ACT[k] + 1].sum()) for i, k in enumerate(CATS) if k in ACT},
                       var_share={k: {c: float(r2["var"][c][i] / sd2[i] ** 2) for c in ("nat", "reg", "seat")} for i, k in enumerate(CATS)})
# ---- lambda scan (2015-size 'shy Tory' miss multiples)
e_scan = FastEngine(m, nat_spec=SCAN_NAT, reg_ns=REG, nat_norm=True)
for lam in (0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5):
    D = m.scenario_shift(lam)
    sh = m.gb_shares(D)
    r = pack(e_scan.run(lam))
    r["gb_shares_centre"] = {k: float(100 * sh[i]) for i, k in enumerate(CATS)}
    res[f"lam{lam}"] = r
    print(f"lam {lam}: centre shares", {k: round(v, 1) for k, v in r["gb_shares_centre"].items()},
          "E seats", {k: round(v["mean"], 1) for k, v in r["summary"].items()}, flush=True)
# ---- oracle: actual 2024 shares as inputs, national factor switched off
POLLS_ACT = {"GB": C.ACTUAL_GB_2024,
             "Scotland": {"SNP": 30.0, "Lab": 35.3, "Con": 12.7, "LD": 9.7, "Grn": 3.8, "Ref": 7.0},
             "Wales": {"Lab": 37.0, "Con": 18.2, "PC": 14.8, "LD": 6.5, "Grn": 4.7, "Ref": 16.9}}
mo = make_model(polls=POLLS_ACT)
eo = FastEngine(mo, nat_spec=[("gh", 1)] * 5, reg_ns=REG)
res["oracle_no_nat"] = pack(eo.run(0.0))
eo2 = FastEngine(mo, nat_spec=SCAN_NAT, reg_ns=REG, nat_norm=True)
res["oracle_with_nat"] = pack(eo2.run(0.0))
print("oracle E seats", {k: round(v["mean"], 1) for k, v in res["oracle_no_nat"]["summary"].items()}, flush=True)
# ---- seat-level calibration (lambda=0 and oracle), moderate rule
act = m.df[[f"a24_{k}" for k in ["Con", "Lab", "LD", "Ref", "Grn", "SNP", "PC", "OthMax"]]].values.argmax(1)
for tag, mod in (("lam0", m), ("oracle", mo)):
    ee = ExactEngine(mod, nat_ns=(9, 5, 3, 2, 2) if tag == "lam0" else (1, 1, 1, 1, 1), reg_ns=(3, 2, 2, 2, 2, 2))
    P = ee.seat_marginals(0.0)
    y = np.zeros_like(P); y[np.arange(len(act)), act] = 1
    brier = float(((P - y) ** 2).sum(1).mean())
    logs = float(-np.log(np.maximum(P[np.arange(len(act)), act], 1e-300)).mean())
    pw = P.max(1); hit = (P.argmax(1) == act)
    bins = np.array([0, .5, .6, .7, .8, .9, .95, .99, 1.0001])
    rel = []
    for lo_, hi_ in zip(bins[:-1], bins[1:]):
        s = (pw >= lo_) & (pw < hi_)
        if s.sum():
            rel.append(dict(lo=float(lo_), hi=float(min(hi_, 1)), n=int(s.sum()), mean_p=float(pw[s].mean()), hit=float(hit[s].mean())))
    # all (seat, party) pairs reliability on predicted probability
    pr = P.ravel(); yr = y.ravel()
    b2 = np.array([0, .01, .05, .2, .4, .6, .8, .95, .99, 1.0001]); rel2 = []
    for lo_, hi_ in zip(b2[:-1], b2[1:]):
        s = (pr >= lo_) & (pr < hi_)
        if s.sum():
            rel2.append(dict(lo=float(lo_), hi=float(min(hi_, 1)), n=int(s.sum()), mean_p=float(pr[s].mean()), freq=float(yr[s].mean())))
    upsets = [dict(seat=m.df.name.iloc[i], region=m.df.region.iloc[i], fav=CATS[int(P[i].argmax())], p_fav=float(pw[i]),
                   winner=CATS[int(act[i])], p_winner=float(P[i, act[i]])) for i in np.where(~hit)[0]]
    upsets.sort(key=lambda d: d["p_winner"])
    res[f"calib_{tag}"] = dict(brier=brier, logscore=logs, hit_rate=float(hit.mean()), expected_hits=float(pw.sum()),
                               reliability_fav=rel, reliability_all=rel2, n_upsets=len(upsets), upsets=upsets)
    print(tag, "Brier %.4f log %.4f hit %.3f expected hits %.1f actual %d" % (brier, logs, hit.mean(), pw.sum(), hit.sum()), flush=True)
json.dump(res, open(f"{OUT}/backtest2_v0_1.json", "w"))
print("done")
