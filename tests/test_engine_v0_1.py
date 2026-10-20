"""Engine tests: (1) Poisson-binomial recursion vs enumeration; (2) a tiny synthetic model where the
engine's whole pipeline (regional mixing + convolution + national weighting) is checked against brute-force
enumeration of every seat outcome at every quadrature node; (3) quadrature rules integrate Gaussian
moments exactly; (4) seat-probability mixture matches Gumbel-max simulation."""
import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import sys, itertools, numpy as np
sys.path.insert(0, _ROOT)
from ukseat.core import exact_chunk
from ukseat.quad import tensor_gh, gh1
from ukseat.engine import outer_rule, trapezoid_normal

rng = np.random.default_rng(7)
ok = True

def check(name, cond, detail=""):
    global ok
    ok &= bool(cond)
    print(("PASS " if cond else "FAIL ") + name + (" " + detail if detail else ""))

# ---------- (2) tiny synthetic model: R=2 regions, 3+4 seats, K=3 parties
K = 3; n = 7; N = n
roff = np.array([0, 3, 7], np.int64)
pres = np.ones((n, K), bool); pres[0, 2] = False          # one party absent in seat 0
base = np.log(rng.dirichlet(np.ones(K), size=n)); base[~pres] = -1e9
pi, t1, t2 = 0.2, 0.11, 0.25
Zr, Wr = tensor_gh((3, 2))
Cr = rng.normal(0, 0.1, size=(2, K, 2))
rs = np.stack([Zr @ Cr[r].T for r in range(2)]); rw = np.stack([Wr, Wr]); rc = np.array([len(Wr)] * 2, np.int64)
Zn, Wn = outer_rule([("gh", 5), ("trap", 0.5, 4.0)])
An = rng.normal(0, 0.2, size=(K, 2))
Dn = Zn @ An.T
d, mu, vs, vr, tv = exact_chunk(Dn, base, pres, roff, rs, rw, rc, pi, t1, t2, 0.0, N)
dist = np.tensordot(Wn, d, 1)

def probs(eta, pr):
    out = np.zeros(K)
    for w, tau in ((1 - pi, t1), (pi, t2)):
        x = np.where(pr, eta / tau, -np.inf); e = np.exp(x - x.max()); out += w * e / e.sum()
    return out

brute = np.zeros((K, N + 1))
for a in range(len(Wn)):
    # regions independent given national node: enumerate each region's outcomes over its nodes
    reg_d = []
    for r in range(2):
        idx = range(roff[r], roff[r + 1])
        rd = np.zeros((K, len(idx) + 1))
        for b in range(rc[r]):
            P = [probs(base[i] + Dn[a] + rs[r, b], pres[i]) for i in idx]
            for outc in itertools.product(range(K), repeat=len(P)):
                pr = np.prod([P[j][outc[j]] for j in range(len(P))])
                for k in range(K):
                    rd[k, sum(1 for o in outc if o == k)] += rw[r, b] * pr
        reg_d.append(rd)
    for k in range(K):
        brute[k] += Wn[a] * np.convolve(reg_d[0][k], reg_d[1][k])[: N + 1]
err = np.abs(brute - dist).max()
check("tiny model: engine == brute-force enumeration", err < 1e-14, f"max abs diff {err:.2e}")
check("tiny model: distributions sum to 1", np.abs(dist.sum(1) - 1).max() < 1e-13)
means = dist @ np.arange(N + 1)
check("tiny model: E[S] from dist == accumulated moment", np.abs(means - Wn @ mu).max() < 1e-12)

# ---------- truncation bound: eps>0 changes the distribution by less than the reported TV bound
d2, _, _, _, tv2 = exact_chunk(Dn, base, pres, roff, rs, rw, rc, pi, t1, t2, 1e-3, N)
dist2 = np.tensordot(Wn, d2, 1)
tvdiff = 0.5 * np.abs(dist2 - dist).sum(1).max()
check("truncation: TV change <= reported bound", tvdiff <= (Wn @ tv2) + 1e-15, f"TV {tvdiff:.2e} <= bound {(Wn @ tv2):.2e}")

# ---------- (3) quadrature rules: E[z^2]=1, E[z^4]=3, E[z^6]=15 per dimension
for n_ in (3, 5, 9):
    x, w = gh1(n_)
    exact_deg = 2 * n_ - 1
    m2, m4, m6 = w @ x ** 2, w @ x ** 4, w @ x ** 6
    check(f"GH{n_} moments", abs(m2 - 1) < 1e-12 and (exact_deg < 4 or abs(m4 - 3) < 1e-11) and (exact_deg < 6 or abs(m6 - 15) < 1e-10))
x, w = trapezoid_normal(0.25, 8.0)
check("trapezoid h=0.25 moments", abs(w @ x ** 2 - 1) < 1e-12 and abs(w @ x ** 4 - 3) < 1e-11)

# ---------- (4) closed-form mixture-softmax vs Gumbel-max simulation (generative model)
eta = np.log(np.array([0.40, 0.35, 0.15, 0.10]))
pr = np.ones(4, bool)
p = np.zeros(4)
for wgt, tau in ((1 - pi, t1), (pi, t2)):
    e = np.exp(eta / tau - (eta / tau).max()); p += wgt * e / e.sum()
M = 2_000_000
tau = np.where(rng.random(M) < pi, t2, t1)
g = -np.log(-np.log(rng.random((M, 4))))
win = np.argmax(eta[None] / tau[:, None] + g, 1)
f = np.bincount(win, minlength=4) / M
z = (f - p) / np.sqrt(p * (1 - p) / M + 1e-300)
check("mixture-softmax == Gumbel-max frequencies", np.abs(z).max() < 4.5, "z=" + np.array2string(z, precision=2))

print("ALL PASS" if ok else "SOME FAILED")
import sys as _sys
_sys.exit(0 if ok else 1)   # R1 review 3.2: nonzero exit on failure
