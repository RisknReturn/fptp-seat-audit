"""Independent generative Monte Carlo for the v0.2 model (numpy; samples z, mu, regional w, seat noise via Gumbel-max)
compared with the QMC-outer / exact-inner engine for the poll-based ex-ante configuration."""
import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import sys, json, numpy as np
sys.path.insert(0, _ROOT)
from ukseat.model2 import SeatModel2
from ukseat.setup2024 import reg_cfg
from ukseat import config_2024 as C
from ukseat.engine2 import QEngine2
from ukseat.report import events, events_mc, summarise_dist
from ukseat.model import CATS
R = _ROOT
m = SeatModel2(f"{R}/data/gb2024_seats_v0_1.csv", f"{R}/data/calib_v0_1.json", C.POLLS_2024, C.SIGMA_NAT, C.CORR_NAT,
               reg_cfg(), e_ref=C.E_REF_2015, cov_csv=f"{R}/data/gb2024_v02_covariates.csv", kappa=4.0)
K = len(CATS); pres = np.isfinite(m.logb[:, :K])
base = np.where(pres, m.logb[:, :K], 0.0) + m.mean_delta()[:, :K]
rng = np.random.default_rng(20261019)
NS, CH = 200_000, 400
S = np.zeros((NS, K), np.int32)
Lk = m.L[:, :K, :]
for s0 in range(0, NS, CH):
    n = min(CH, NS - s0)
    z = rng.standard_normal((n, m.A.shape[1])); dnat = z @ m.A[:K].T                       # (n,K)
    mu = m.mu_mean + m.mu_sd * rng.standard_normal((n, len(m.mu_mean)))                    # (n,J)
    eta = base[None] + dnat[:, None, :] + np.einsum("ikj,sj->sik", Lk, mu)
    for j, Cj in enumerate(m.C):
        idx = m.reg == j
        w = rng.standard_normal((n, Cj.shape[1]))
        eta[:, idx, :] += (w @ Cj[:K].T)[:, None, :]
    tau = np.where(rng.random((n, m.n)) < m.pi, m.tau2, m.tau1)
    g = -np.log(-np.log(rng.random((n, m.n, K))))
    v = np.where(pres[None], eta / tau[:, :, None] + g, -np.inf)
    win = v.argmax(2)
    for k in range(K):
        S[s0:s0 + n, k] = (win == k).sum(1)
evm = events_mc(S)
q = QEngine2(m).run(m_log2=11, reps=4)
eve = events(q["dist"]); evr = [events(x) for x in q["dist_reps"]]
out = {}
print(f"{'event':28s} {'QMC+exact':>10s} {'(se)':>8s} {'MC 2e5':>9s} {'z':>6s}")
for k in eve:
    se_q = np.std([e[k] for e in evr], ddof=1) / np.sqrt(len(evr))
    se_m = np.sqrt(max(evm[k] * (1 - evm[k]), 1e-12) / NS)
    zz = (eve[k] - evm[k]) / np.sqrt(se_q ** 2 + se_m ** 2)
    out[k] = dict(qmc=eve[k], qmc_se=se_q, mc=evm[k], z=zz)
    print(f"{k:28s} {eve[k]:10.4f} {se_q:8.4f} {evm[k]:9.4f} {zz:6.1f}")
mq = summarise_dist(q["dist"]); print("means QMC", {k: round(v["mean"], 2) for k, v in mq.items()}); print("means MC ", dict(zip(CATS, S.mean(0).round(2))))
json.dump(out, open(f"{R}/out/v02_mc_check.json", "w"), indent=1)
