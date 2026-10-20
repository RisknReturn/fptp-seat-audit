"""Exact seat-distribution engine and Monte Carlo baseline.

Exact: outer Gauss-Hermite over national factor z; for each region, inner Gauss-Hermite over the
regional factor w_r; given (z, w_r) seats are independent, so each party's regional seat count is
Poisson-binomial (computed by the Andersen-Sidenius-Basu style recursion). Regional distributions
are averaged over w_r, convolved across regions (independent given z), and averaged over z.
No sampling anywhere. Given the factors, the recursion and convolution are algebraically exact; the
implementation adds floating-point error and the optional truncation of near-certain seats (bounded in total
variation). The dominant approximation is the quadrature over the factors (checked by refinement; not a bound).
"""
import time
import numpy as np
from numba import njit, prange
from .model import ALL, CATS, IDX
from .quad import tensor_gh

KC = len(CATS)


@njit(cache=True)
def _pb(p):
    n = p.shape[0]
    d = np.zeros(n + 1)
    d[0] = 1.0
    for i in range(n):
        q = p[i]
        for s in range(i + 1, 0, -1):
            d[s] = d[s] * (1.0 - q) + d[s - 1] * q
        d[0] *= 1.0 - q
    return d


@njit(parallel=True, cache=True)
def region_dist(P, wts, eps, nmax):
    """P: (m, n, K) win probs at m regional nodes; wts: (m,).
    Returns dist (K, nmax+1) averaged over nodes, and per-party TV truncation bound (K,).
    Seats with p<=eps at every node are dropped; p>=1-eps at every node are fixed wins."""
    m, n, K = P.shape
    out = np.zeros((K, nmax + 1))
    tv = np.zeros(K)
    for k in prange(K):
        # classify seats for party k
        act = np.zeros(n, np.int64)
        na = 0
        fixed = 0
        b = 0.0
        for i in range(n):
            lo = 1.0; hi = 0.0
            for a in range(m):
                v = P[a, i, k]
                if v < lo: lo = v
                if v > hi: hi = v
            if hi <= eps:
                b += hi
            elif lo >= 1.0 - eps:
                fixed += 1
                b += 1.0 - lo
            else:
                act[na] = i
                na += 1
        tv[k] = b
        pa = np.empty(na)
        for a in range(m):
            for j in range(na):
                pa[j] = P[a, act[j], k]
            d = _pb(pa)
            for s in range(na + 1):
                out[k, s + fixed] += wts[a] * d[s]
    return out, tv


def _rotate_by_impact(model, G):
    """Orthogonal rotation of a standard-normal factor so dim 0 has the largest seat impact."""
    M = G.T @ G
    ev, Q = np.linalg.eigh(M)
    o = np.argsort(ev)[::-1]
    return Q[:, o], ev[o]


class ExactEngine:
    def __init__(self, model, nat_ns=(9, 5, 3, 3, 2), reg_ns=(5, 3, 2, 2, 2, 2), eps=1e-10, rotate=True):
        self.m = model
        self.nat_ns, self.eps = tuple(nat_ns), eps
        self.reg_ns = reg_ns if isinstance(reg_ns, (str, dict)) else tuple(reg_ns)
        self.N = model.n
        self.present = np.isfinite(model.logb[:, :KC])
        self.base = np.where(self.present, model.logb[:, :KC], -np.inf) + model.mean_delta()[:, :KC]
        self.rid = [np.where(model.reg == j)[0] for j in range(len(model.regions))]
        # rotations by seat impact at the centre (allocation of nodes only; integral unchanged)
        self.Qn = np.eye(model.A.shape[1])
        self.Qr = [np.eye(C.shape[1]) for C in model.C]
        if rotate:
            self._set_rotations()

    def _impact(self, Delta, cols):
        """d E[seats_k]/d factor for loading matrix cols (K_all x d), numerically."""
        h = 1e-4
        def es(D):
            eta = self.base + D[None, :KC]
            return self.m.win_probs(eta, self.present).sum(0)
        f0 = es(Delta)
        G = np.zeros((KC, cols.shape[1]))
        for j in range(cols.shape[1]):
            G[:, j] = (es(Delta + h * cols[:, j]) - f0) / h
        return G

    def _set_rotations(self):
        D0 = np.zeros(len(ALL))
        Gn = self._impact(D0, self.m.A)
        if getattr(self, "nat_norm", False):
            # party-normalised impact (as for regions): small parties' drivers are not ranked last
            P0 = self.m.win_probs(self.base + D0[None, :KC], self.present)
            Gn = Gn / np.sqrt(np.maximum((P0 * (1 - P0)).sum(0), 0.25))[:, None]
        self.Qn, self.nat_ev = _rotate_by_impact(self.m, Gn)
        self.Qr = []; self.reg_ev = []
        for j, C in enumerate(self.m.C):
            idx = self.rid[j]
            def es(D, idx=idx):
                eta = self.base[idx] + D[None, :KC]
                return self.m.win_probs(eta, self.present[idx]).sum(0)
            h = 1e-4; f0 = es(D0)
            G = np.zeros((KC, C.shape[1]))
            for c in range(C.shape[1]):
                G[:, c] = (es(D0 + h * C[:, c]) - f0) / h
            # party-normalised impact: each party's seat response scaled by its seat sd in the region
            P0 = self.m.win_probs(self.base[idx] + D0[None, :KC], self.present[idx])
            sk = np.sqrt(np.maximum((P0 * (1 - P0)).sum(0), 0.25))
            Q, ev = _rotate_by_impact(self.m, G / sk[:, None])
            self.Qr.append(Q); self.reg_ev.append(ev)

    def reg_ns_for(self, j):
        C = self.m.C[j]
        if self.reg_ns == "auto" or isinstance(self.reg_ns, dict):
            cfg = {"hi": 0.5, "mid": 0.1, "lo": 0.01} if self.reg_ns == "auto" else self.reg_ns
            ev = self.reg_ev[j]; sh = ev / ev.sum()
            return tuple(5 if x > cfg["hi"] else 3 if x > cfg["mid"] else 2 if x > cfg["lo"] else 1 for x in sh)
        return tuple(self.reg_ns[: C.shape[1]])

    def run(self, lam=0.0, keep_conditional=False):
        t0 = time.time()
        m = self.m
        Dl = m.scenario_shift(lam)
        Zn, Wn = tensor_gh(self.nat_ns[: m.A.shape[1]])
        Zn = Zn @ self.Qn.T
        reg_rules = []
        for j, C in enumerate(m.C):
            X, W = tensor_gh(self.reg_ns_for(j))
            reg_rules.append(((X @ self.Qr[j].T) @ C[:KC].T, W))   # (m_r, KC) delta shifts
        N = self.N
        dist = np.zeros((KC, N + 1))
        mom = np.zeros((3, KC))            # E[S], E[Var(S|z)], E[Var(S|z,w)] pieces
        cond_mean = np.zeros((len(Wn), KC))
        tvmax = 0.0
        for a in range(len(Wn)):
            Dn = Dl + m.A @ Zn[a]
            tot = np.zeros((KC, 1)); tot[:, 0] = 1.0
            var_seat = np.zeros(KC); var_reg = np.zeros(KC); mean_z = np.zeros(KC)
            tv = np.zeros(KC)
            for j, idx in enumerate(self.rid):
                shifts, Wr = reg_rules[j]
                eta = self.base[idx][None, :, :] + Dn[None, None, :KC] + shifts[:, None, :]
                P = m.win_probs(eta, self.present[idx][None])
                d, tvr = region_dist(P, Wr, self.eps, len(idx))
                tv += tvr
                mu_w = P.sum(1)                                   # (m_r, KC)
                mu = Wr @ mu_w
                var_seat += Wr @ (P * (1 - P)).sum(1)
                var_reg += Wr @ (mu_w ** 2) - mu ** 2
                mean_z += mu
                tot = np.array([np.convolve(tot[k], d[k]) for k in range(KC)])
            dist += Wn[a] * tot[:, : N + 1]
            mom[0] += Wn[a] * mean_z
            mom[1] += Wn[a] * var_seat
            mom[2] += Wn[a] * var_reg
            cond_mean[a] = mean_z
            tvmax = max(tvmax, tv.max())
        var_nat = Wn @ (cond_mean ** 2) - mom[0] ** 2
        res = dict(dist=dist, mean=mom[0], var_seat=mom[1], var_reg=mom[2], var_nat=var_nat,
                   tv_bound=tvmax, nodes_outer=len(Wn), seconds=time.time() - t0, lam=lam)
        if keep_conditional:
            res["cond_mean"], res["Zn"], res["Wn"] = cond_mean, Zn, Wn
        return res

    def seat_marginals(self, lam=0.0):
        """E[p_ik] over factors (for seat-level calibration), using the same rules."""
        m = self.m
        Dl = m.scenario_shift(lam)
        Zn, Wn = tensor_gh(self.nat_ns[: m.A.shape[1]]); Zn = Zn @ self.Qn.T
        out = np.zeros((self.N, KC))
        for j, C in enumerate(m.C):
            idx = self.rid[j]
            X, W = tensor_gh(self.reg_ns_for(j))
            sh = (X @ self.Qr[j].T) @ C[:KC].T
            for a in range(len(Wn)):
                Dn = Dl + m.A @ Zn[a]
                eta = self.base[idx][None] + Dn[None, None, :KC] + sh[:, None, :]
                out[idx] += Wn[a] * np.tensordot(W, m.win_probs(eta, self.present[idx][None]), 1)
        return out


# ------------------------------------------------------------------ Monte Carlo baseline
@njit(parallel=True, cache=True)
def _mc_kernel(base, present, reg, A, Cs, Dl, pi, t1, t2, nsim, seed, out):
    """Simulates the generative model directly: Gaussian factors + Gumbel-max seat noise.
    base: (n,KC) eta at centre; Cs: (R, KC, dmax) regional loadings; out: (nsim, KC) seat counts."""
    n, K = base.shape
    R = Cs.shape[0]
    dn = A.shape[1]
    dr = Cs.shape[2]
    for s in prange(nsim):
        np.random.seed(seed + s)
        z = np.random.standard_normal(dn)
        dnat = Dl.copy()
        for k in range(K):
            for c in range(dn):
                dnat[k] += A[k, c] * z[c]
        dreg = np.zeros((R, K))
        for r in range(R):
            w = np.random.standard_normal(dr)
            for k in range(K):
                acc = 0.0
                for c in range(dr):
                    acc += Cs[r, k, c] * w[c]
                dreg[r, k] = acc
        for i in range(n):
            tau = t2 if np.random.random() < pi else t1
            best = -1e300; bk = -1
            for k in range(K):
                if present[i, k]:
                    u = np.random.random()
                    g = -np.log(-np.log(u))
                    v = (base[i, k] + dnat[k] + dreg[reg[i], k]) / tau + g
                    if v > best:
                        best = v; bk = k
            out[s, bk] += 1


def monte_carlo(model, lam, nsim, seed=12345):
    t0 = time.time()
    present = np.isfinite(model.logb[:, :KC])
    base = np.where(present, model.logb[:, :KC], 0.0) + model.mean_delta()[:, :KC]
    dmax = max(C.shape[1] for C in model.C)
    Cs = np.zeros((len(model.C), KC, dmax))
    for j, C in enumerate(model.C):
        Cs[j, :, : C.shape[1]] = C[:KC]
    out = np.zeros((nsim, KC), np.int32)
    Dl = model.scenario_shift(lam)[:KC].copy()
    _mc_kernel(base, present, model.reg.astype(np.int64), model.A[:KC].copy(), Cs, Dl,
               model.pi, model.tau1, model.tau2, nsim, seed, out)
    return out, time.time() - t0


# ------------------------------------------------------------------ fast engine (numba core)
from .core import exact_chunk
from .quad import gh1


def trapezoid_normal(h, L):
    x = np.arange(-L, L + 1e-12, h)
    w = np.exp(-0.5 * x ** 2); w /= w.sum()
    return x, w


def outer_rule(spec):
    """spec: list per rotated national dim: ('gh', n) or ('trap', h, L). Returns nodes (M,d), weights."""
    import itertools
    rules = []
    for s in spec:
        if s[0] == "gh":
            rules.append(gh1(s[1]))
        else:
            rules.append(trapezoid_normal(s[1], s[2]))
    X = np.array(list(itertools.product(*[r[0] for r in rules])))
    W = np.array([np.prod(c) for c in itertools.product(*[r[1] for r in rules])])
    return X.reshape(len(W), len(spec)), W


class FastEngine(ExactEngine):
    """Same integral as ExactEngine, compiled and parallel over national nodes."""

    def __init__(self, model, nat_spec, reg_ns=(5, 3, 2, 2, 2, 2), eps=1e-12, chunk=256, rotate=True, nat_norm=False):
        self.nat_norm = nat_norm
        super().__init__(model, nat_ns=(1,) * model.A.shape[1], reg_ns=reg_ns, eps=eps, rotate=rotate)
        self.nat_spec, self.chunk = nat_spec, chunk
        order = np.argsort(model.reg, kind="stable")
        self.order = order
        self.base_s = np.where(self.present, self.base, -1e9)[order]
        self.pres_s = self.present[order]
        counts = np.bincount(model.reg, minlength=len(model.regions))
        self.roff = np.concatenate([[0], np.cumsum(counts)]).astype(np.int64)

    def _reg_arrays(self):
        rules = []
        for j, C in enumerate(self.m.C):
            X, W = tensor_gh(self.reg_ns_for(j))
            rules.append(((X @ self.Qr[j].T) @ C[:KC].T, W))
        Mmax = max(len(r[1]) for r in rules)
        R = len(rules)
        rs = np.zeros((R, Mmax, KC)); rw = np.zeros((R, Mmax)); rc = np.zeros(R, np.int64)
        for j, (S, W) in enumerate(rules):
            rs[j, : len(W)] = S; rw[j, : len(W)] = W; rc[j] = len(W)
        return rs, rw, rc

    def run(self, lam=0.0, keep_conditional=False):
        t0 = time.time()
        m = self.m
        Dl = m.scenario_shift(lam)
        Z, W = outer_rule(self.nat_spec)
        Dn = (Dl[None, :] + (Z @ self.Qn.T) @ m.A.T)[:, :KC].copy()
        rs, rw, rc = self._reg_arrays()
        N = self.N
        dist = np.zeros((KC, N + 1)); mean = np.zeros(KC); vseat = np.zeros(KC); vreg = np.zeros(KC)
        cm = np.zeros((len(W), KC)); tvb = 0.0
        for s in range(0, len(W), self.chunk):
            sl = slice(s, s + self.chunk)
            d, mu, vs, vr, tv = exact_chunk(Dn[sl], self.base_s, self.pres_s, self.roff, rs, rw, rc,
                                            m.pi, m.tau1, m.tau2, self.eps, N)
            w = W[sl]
            dist += np.tensordot(w, d, 1); mean += w @ mu; vseat += w @ vs; vreg += w @ vr
            cm[sl] = mu; tvb = max(tvb, tv.max())
        var_nat = W @ (cm ** 2) - mean ** 2
        res = dict(dist=dist, mean=mean, var_seat=vseat, var_reg=vreg, var_nat=var_nat, tv_bound=tvb,
                   nodes_outer=len(W), nodes_reg=[int(x) for x in rc], seconds=time.time() - t0, lam=lam)
        if keep_conditional:
            res.update(cond_mean=cm, Z=Z, W=W)
        return res
