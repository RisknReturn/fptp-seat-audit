"""KAPPA sensitivity (Oth loading on mu_L) for the oracle run with realised mu; seat-level check of Oth winners."""
import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import sys, json, numpy as np
sys.path.insert(0, _ROOT)
sys.argv = [sys.argv[0], "4"]
src = open(f"{_ROOT}/scripts/v02_experiments.py").read()
exec(src[:src.index("out = {\"kappa\"")])  # helpers incl. Raked
import ukseat.model2 as M2
res = {}
for kap in (0.0, 2.0, 3.0, 4.0, 6.0):
    KAPPA = kap
    mo = SeatModel2(f"{ROOT}/data/gb2024_seats_v0_1.csv", f"{ROOT}/data/calib_v0_1.json", POLLS_ACT, C.SIGMA_NAT, C.CORR_NAT,
                    reg_cfg(), e_ref=C.E_REF_2015, cov_csv=f"{ROOT}/data/gb2024_v02_covariates.csv", kappa=kap)
    with Raked(mo, MU_REAL):
        r = QEngine2(mo, nat_on=False, mu_fixed=MU_REAL).run(reps=1, m_log2=0)
        s = summarise_dist(r["dist"])
        P = seat_marginals(mo, 0.0, False, MU_REAL, ndraw=400)
    oth = P[:, CATS.index("Oth")]
    top = np.argsort(-oth)[:10]
    act = mo.df[[f"a24_{k}" for k in ["Con", "Lab", "LD", "Ref", "Grn", "SNP", "PC", "OthMax"]]].values.argmax(1)
    res[kap] = dict(Lab=s["Lab"]["mean"], Oth=s["Oth"]["mean"], Oth_q=(s["Oth"]["q05"], s["Oth"]["q95"]), cal=calib(mo, P),
                    top_oth=[(str(mo.df["name"].iloc[i]), round(float(oth[i]), 2), bool(act[i] == 7)) for i in top])
    print(f"kappa {kap}: Lab {s['Lab']['mean']:.1f} Oth {s['Oth']['mean']:.1f} ({s['Oth']['q05']}-{s['Oth']['q95']}) brier {res[kap]['cal']['brier']:.4f}", flush=True)
    print("   top P(Oth):", [(n[:22], p, w) for n, p, w in res[kap]["top_oth"][:8]], flush=True)
json.dump({str(k): v for k, v in res.items()}, open(f"{ROOT}/out/kappa_sens_v02.json", "w"), indent=1, default=str)
