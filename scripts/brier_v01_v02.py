"""Seat-level calibration of v0.1 (no nominations, no mu) vs v0.2, oracle and poll inputs, same scoring code."""
import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import sys, json, numpy as np
sys.path.insert(0, _ROOT)
sys.argv = [sys.argv[0], "4"]
exec(open(f"{_ROOT}/scripts/v02_experiments.py").read().split("out = {")[0])
from ukseat.model import SeatModel
def mk1(polls):
    return SeatModel(f"{ROOT}/data/gb2024_seats_v0_1.csv", f"{ROOT}/data/calib_v0_1.json", polls, C.SIGMA_NAT, C.CORR_NAT,
                     reg_cfg(), e_ref=C.E_REF_2015)
def marg1(model, nat_on, ndraw=1500, seed=3):
    rng = np.random.default_rng(seed); K = len(CATS); pres = np.isfinite(model.logb[:, :K])
    base = np.where(pres, model.logb[:, :K], 0.0) + model.mean_delta()[:, :K]
    P = np.zeros((model.n, K))
    for _ in range(ndraw):
        eta = base.copy()
        if nat_on: eta += (model.A @ rng.standard_normal(model.A.shape[1]))[:K]
        for j, C_ in enumerate(model.C):
            idx = model.reg == j; eta[idx] += (C_ @ rng.standard_normal(C_.shape[1]))[:K]
        P += model.win_probs(eta, pres)
    return P / ndraw
out = {}
m1o, m1p = mk1(POLLS_ACT), mk1(C.POLLS_2024)
out["v01_oracle"] = calib(m1o, marg1(m1o, False)); print("v01 oracle", out["v01_oracle"], flush=True)
out["v01_polls"] = calib(m1p, marg1(m1p, True)); print("v01 polls ", out["v01_polls"], flush=True)
json.dump(out, open(f"{ROOT}/out/brier_v01_v02.json", "w"), indent=1)
