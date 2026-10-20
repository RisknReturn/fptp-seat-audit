"""v0.2 engine tests (review of 2026-10-21):
(1) the v0.2 core exact_chunk_L (seat-type factor term L @ mu) against brute-force enumeration on a tiny model;
(2) QEngine2.seat_marginals() is consistent with QEngine2.run(): for the same arguments, the sum over seats of the
    seat-level win probabilities equals the mean of the seat distribution, party by party and scrambling by
    scrambling, on the bundled final-poll model (national factors on/off, mu fixed at the prior mean or the realised
    2024 values, normal and Student-t mu, the constructor default and the audited regional rule);
(3) the same consistency for FastEngine (v0.1 model, outer rule nat_spec).
The inherited ExactEngine method failed (2): it omitted L @ mu and used one national node (about 4.8 seats for Labour
with mu at the prior mean, 33 seats with the realised mu)."""
import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import sys, json, itertools, numpy as np
sys.path.insert(0, _ROOT)
from ukseat.core import exact_chunk_L
from ukseat.quad import tensor_gh
from ukseat.engine import FastEngine, KC
from ukseat.engine2 import QEngine2, REG_NS_AUDITED
from ukseat.model2 import SeatModel2, FACTORS
from ukseat.setup2024 import reg_cfg, make_model
from ukseat import config_2024 as C

rng = np.random.default_rng(11)
ok = True

def check(name, cond, detail=""):
    global ok
    ok &= bool(cond)
    print(("PASS " if cond else "FAIL ") + name + (" " + detail if detail else ""), flush=True)

# ---------- (1) exact_chunk_L vs enumeration: R = 2 regions, 3 + 3 seats, K = 3 parties, J = 2 seat-type factors
K = 3; n = 6; N = n; J = 2
roff = np.array([0, 3, 6], np.int64)
pres = np.ones((n, K), bool); pres[1, 2] = False
base = np.log(rng.dirichlet(np.ones(K), size=n)); base[~pres] = -1e9
Lseat = rng.normal(0, 0.15, size=(n, K, J)); Lseat[~pres] = 0.0
pi, t1, t2 = 0.2, 0.11, 0.25
Zr, Wr = tensor_gh((3, 2))
Cr = rng.normal(0, 0.1, size=(2, K, 2))
rs = np.stack([Zr @ Cr[r].T for r in range(2)]); rw = np.stack([Wr, Wr]); rc = np.array([len(Wr)] * 2, np.int64)
A_ = 5
Dn = rng.normal(0, 0.1, size=(A_, K)); Mu = rng.normal(0, 1.0, size=(A_, J))
d, _, _, _, _ = exact_chunk_L(Dn, base, pres, Lseat, Mu, roff, rs, rw, rc, pi, t1, t2, 0.0, N)

def probs(eta, pr):
    out = np.zeros(K)
    for w, tau in ((1 - pi, t1), (pi, t2)):
        x = np.where(pr, eta / tau, -np.inf); e = np.exp(x - x.max()); out += w * e / e.sum()
    return out

err = 0.0
for a in range(A_):
    reg_d = []
    for r in range(2):
        idx = range(roff[r], roff[r + 1])
        rd = np.zeros((K, len(idx) + 1))
        for b in range(rc[r]):
            P = [probs(base[i] + Dn[a] + rs[r, b] + Lseat[i] @ Mu[a], pres[i]) for i in idx]
            for outc in itertools.product(range(K), repeat=len(P)):
                pr = np.prod([P[j][outc[j]] for j in range(len(P))])
                for k in range(K):
                    rd[k, sum(1 for o in outc if o == k)] += rw[r, b] * pr
        reg_d.append(rd)
    brute = np.array([np.convolve(reg_d[0][k], reg_d[1][k])[: N + 1] for k in range(K)])
    err = max(err, np.abs(brute - d[a]).max())
check("exact_chunk_L (with L @ mu) == brute-force enumeration", err < 1e-14, f"max abs diff {err:.2e}")

# ---------- (2) QEngine2: sum_i seat_marginals == E[S] from run(), same arguments
mp = SeatModel2(f"{_ROOT}/data/gb2024_seats_v0_1.csv", f"{_ROOT}/data/calib_v0_1.json", C.POLLS_2024, C.SIGMA_NAT,
                C.CORR_NAT, reg_cfg(), e_ref=C.E_REF_2015, cov_csv=f"{_ROOT}/data/gb2024_v02_covariates.csv",
                kappa=4.0, neutral=True)
REAL = json.load(open(f"{_ROOT}/data/realised_mu_2024.json"))
MU_REAL = np.array([REAL[f]["realised_2024"] for f in FACTORS])
TOL = 1e-7          # seats; run() drops near-certain seats (eps = 1e-12, at most 650 x 1e-12 per party and point)
s = np.arange(mp.n + 1)
cases = [
    ("no national factor, mu fixed at prior mean, regional rule (1,...,1)",
     dict(nat_on=False, mu_fixed=mp.mu_mean, reg_ns=(1,) * 6), dict(m_log2=0, reps=1)),
    ("no national factor, mu fixed at realised 2024 values, regional rule (1,...,1)",
     dict(nat_on=False, mu_fixed=MU_REAL, reg_ns=(1,) * 6), dict(m_log2=0, reps=1)),
    ("all factors, normal mu, regional rule (2,2,1,1,1,1), 2^4 x 2",
     dict(reg_ns=(2, 2, 1, 1, 1, 1)), dict(m_log2=4, reps=2, seed=5)),
    ("all factors, Student-t mu (df 2, scale sqrt(4/3)), lambda 0.5, 2^3 x 2",
     dict(reg_ns=(2, 1, 1, 1, 1, 1)), dict(lam=0.5, m_log2=3, reps=2, seed=9, mu_df=2, mu_scale=np.sqrt(4 / 3))),
    ("constructor default regional rule, mu fixed at prior mean, one outer point",
     dict(nat_on=False, mu_fixed=mp.mu_mean), dict(m_log2=0, reps=1)),
    ("audited regional rule (3,3,3,9,5,3), mu fixed at prior mean, one outer point",
     dict(nat_on=False, mu_fixed=mp.mu_mean, reg_ns=REG_NS_AUDITED), dict(m_log2=0, reps=1)),
]
for name, ekw, rkw in cases:
    e = QEngine2(mp, **ekw)
    r = e.run(**rkw)
    P = e.seat_marginals(return_reps=True, **rkw)
    es_run = r["dist_reps"] @ s                       # (reps, KC)
    es_mar = P.sum(1)                                 # (reps, KC)
    dev = np.abs(es_run - es_mar).max()
    check(f"QEngine2 marginals vs run: {name}", dev < TOL and np.all((P >= 0) & (P <= 1 + 1e-12)),
          f"max |diff| {dev:.1e} seats (Lab {es_run.mean(0)[1]:.4f})")
# outer points of seat_marginals are those of run(): a different seed must move the expected seats
e = QEngine2(mp, reg_ns=(1,) * 6)
a = e.seat_marginals(m_log2=2, reps=1, seed=0).sum(0); b = e.run(m_log2=2, reps=1, seed=1)["dist"] @ s
check("QEngine2 marginals depend on the outer seed (control)", np.abs(a - b).max() > 1e-3, f"max |diff| {np.abs(a - b).max():.2f}")

# ---------- (3) FastEngine: sum_i seat_marginals == run() mean (v0.1 model)
m1 = make_model()
d1 = m1.A.shape[1]
fe = FastEngine(m1, nat_spec=[("gh", 3), ("gh", 2)] + [("gh", 1)] * (d1 - 2), reg_ns=(2, 1, 1, 1, 1, 1))
r1 = fe.run(0.0)
P1 = fe.seat_marginals(0.0)
dev = max(np.abs(P1.sum(0) - r1["mean"]).max(), np.abs(P1.sum(0) - r1["dist"] @ np.arange(m1.n + 1)).max())
check("FastEngine marginals vs run (outer rule nat_spec)", dev < TOL, f"max |diff| {dev:.1e} seats")

print("ALL PASS" if ok else "SOME FAILED")
sys.exit(0 if ok else 1)
