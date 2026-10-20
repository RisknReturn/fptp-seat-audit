"""Numba core of the exact engine: parallel over national quadrature nodes.

For each national node a:
  for each region r:
     for each regional node b:  p_ik = mixture-softmax(eta_ik / tau) for seats i in r
        per party k: Poisson-binomial recursion over seats with eps < p < 1-eps
        (p >= 1-eps -> fixed win; p <= eps -> dropped; dropped mass is added to the TV bound)
     regional distribution = sum_b w_b * dist_b
  national-node distribution = convolution over regions (independent given z)
Returns per-node distributions (chunk, K, N+1) and moments; the caller applies node weights.
"""
import numpy as np
from numba import njit, prange


@njit(parallel=True, cache=True)
def exact_chunk(Dn, base, pres, roff, rshift, rw, rcount, pi, t1, t2, eps, N):
    """exp(eta/tau) = exp(base/tau) * exp(Dn/tau) * exp(rshift/tau): exponentials are factorised so the
    inner loop multiplies instead of calling exp (base = log share + block mean shift; magnitudes stay far from overflow)."""
    A, K = Dn.shape
    n = base.shape[0]
    E1 = np.zeros((n, K)); E2 = np.zeros((n, K))
    for i in range(n):
        for k in range(K):
            if pres[i, k]:
                E1[i, k] = np.exp(base[i, k] / t1); E2[i, k] = np.exp(base[i, k] / t2)
    return _exact_chunk(Dn, E1, E2, roff, rshift, rw, rcount, pi, t1, t2, eps, N)


@njit(parallel=True, cache=True)
def _exact_chunk(Dn, E1, E2, roff, rshift, rw, rcount, pi, t1, t2, eps, N):
    """Dn: (A, K) national delta shifts; base: (n, K) centre eta (seats sorted by region);
    pres: (n, K) bool; roff: (R+1,) region offsets; rshift: (R, Mmax, K); rw: (R, Mmax); rcount: (R,).
    Returns dist (A, K, N+1), mean (A, K), var_seat (A, K), var_reg (A, K), tv (A,)."""
    A, K = Dn.shape
    R = roff.shape[0] - 1
    dist = np.zeros((A, K, N + 1))
    mean = np.zeros((A, K)); vseat = np.zeros((A, K)); vreg = np.zeros((A, K)); tv = np.zeros(A)
    for a in prange(A):
        tot = np.zeros((K, N + 1)); tot[:, 0] = 1.0
        g1 = np.empty(K); g2 = np.empty(K); h1 = np.empty(K); h2 = np.empty(K)
        for k in range(K):
            g1[k] = np.exp(Dn[a, k] / t1); g2[k] = np.exp(Dn[a, k] / t2)
        lo = np.zeros(K, np.int64); hi = np.zeros(K, np.int64)     # support of tot[k]
        eta = np.empty(K)
        for r in range(R):
            i0 = roff[r]; i1 = roff[r + 1]; nr = i1 - i0
            rd = np.zeros((K, nr + 1))
            rlo = np.full(K, nr + 1, np.int64); rhi = np.zeros(K, np.int64)
            P = np.empty((nr, K))
            mu_r = np.zeros(K); m2_r = np.zeros(K)
            tvr = 0.0
            for b in range(rcount[r]):
                wb = rw[r, b]
                for k in range(K):
                    h1[k] = g1[k] * np.exp(rshift[r, b, k] / t1); h2[k] = g2[k] * np.exp(rshift[r, b, k] / t2)
                for i in range(nr):
                    s1 = 0.0; s2 = 0.0
                    for k in range(K):
                        s1 += E1[i0 + i, k] * h1[k]; s2 += E2[i0 + i, k] * h2[k]
                    c1 = (1.0 - pi) / s1; c2 = pi / s2
                    for k in range(K):
                        P[i, k] = c1 * E1[i0 + i, k] * h1[k] + c2 * E2[i0 + i, k] * h2[k]
                for k in range(K):
                    # classify seats, then recursion over active ones
                    fixed = 0; na = 0; drop = 0.0; mu = 0.0; vs = 0.0
                    act = np.empty(nr)
                    for i in range(nr):
                        q = P[i, k]
                        mu += q; vs += q * (1.0 - q)
                        if q <= eps:
                            drop += q
                        elif q >= 1.0 - eps:
                            fixed += 1; drop += 1.0 - q
                        else:
                            act[na] = q; na += 1
                    tvr += wb * drop
                    mu_r[k] += wb * mu; m2_r[k] += wb * mu * mu
                    vseat[a, k] += wb * vs
                    d = np.zeros(na + 1); d[0] = 1.0
                    for j in range(na):
                        q = act[j]
                        for s in range(j + 1, 0, -1):
                            d[s] = d[s] * (1.0 - q) + d[s - 1] * q
                        d[0] *= 1.0 - q
                    for s in range(na + 1):
                        rd[k, s + fixed] += wb * d[s]
                    if fixed < rlo[k]:
                        rlo[k] = fixed
                    if fixed + na > rhi[k]:
                        rhi[k] = fixed + na
            tv[a] += tvr
            for k in range(K):
                mean[a, k] += mu_r[k]
                vreg[a, k] += m2_r[k] - mu_r[k] * mu_r[k]
                # convolve tot[k][lo:hi+1] with rd[k][rlo:rhi+1]
                nlo = lo[k] + rlo[k]; nhi = hi[k] + rhi[k]
                new = np.zeros(nhi - nlo + 1)
                for x in range(lo[k], hi[k] + 1):
                    tx = tot[k, x]
                    if tx == 0.0:
                        continue
                    for y in range(rlo[k], rhi[k] + 1):
                        new[x + y - nlo] += tx * rd[k, y]
                for x in range(lo[k], hi[k] + 1):
                    tot[k, x] = 0.0
                for x in range(nlo, nhi + 1):
                    tot[k, x] = new[x - nlo]
                lo[k] = nlo; hi[k] = nhi
        for k in range(K):
            for x in range(N + 1):
                dist[a, k, x] = tot[k, x]
    return dist, mean, vseat, vreg, tv


@njit(parallel=True, cache=True)
def _exact_chunk_L(Dn, E1, E2, Lseat, Mu, roff, rshift, rw, rcount, pi, t1, t2, eps, N):
    """Dn: (A, K) national delta shifts; base: (n, K) centre eta (seats sorted by region);
    pres: (n, K) bool; roff: (R+1,) region offsets; rshift: (R, Mmax, K); rw: (R, Mmax); rcount: (R,).
    Returns dist (A, K, N+1), mean (A, K), var_seat (A, K), var_reg (A, K), tv (A,)."""
    A, K = Dn.shape
    R = roff.shape[0] - 1
    dist = np.zeros((A, K, N + 1))
    mean = np.zeros((A, K)); vseat = np.zeros((A, K)); vreg = np.zeros((A, K)); tv = np.zeros(A)
    for a in prange(A):
        tot = np.zeros((K, N + 1)); tot[:, 0] = 1.0
        g1 = np.empty(K); g2 = np.empty(K); h1 = np.empty(K); h2 = np.empty(K)
        for k in range(K):
            g1[k] = np.exp(Dn[a, k] / t1); g2[k] = np.exp(Dn[a, k] / t2)
        lo = np.zeros(K, np.int64); hi = np.zeros(K, np.int64)     # support of tot[k]
        eta = np.empty(K)
        for r in range(R):
            i0 = roff[r]; i1 = roff[r + 1]; nr = i1 - i0
            rd = np.zeros((K, nr + 1))
            rlo = np.full(K, nr + 1, np.int64); rhi = np.zeros(K, np.int64)
            P = np.empty((nr, K))
            X1 = np.empty((nr, K)); X2 = np.empty((nr, K))
            J = Mu.shape[1]
            for i in range(nr):
                for k in range(K):
                    ex = 0.0
                    for j in range(J):
                        ex += Lseat[i0 + i, k, j] * Mu[a, j]
                    X1[i, k] = np.exp(ex / t1); X2[i, k] = np.exp(ex / t2)
            mu_r = np.zeros(K); m2_r = np.zeros(K)
            tvr = 0.0
            for b in range(rcount[r]):
                wb = rw[r, b]
                for k in range(K):
                    h1[k] = g1[k] * np.exp(rshift[r, b, k] / t1); h2[k] = g2[k] * np.exp(rshift[r, b, k] / t2)
                for i in range(nr):
                    s1 = 0.0; s2 = 0.0
                    for k in range(K):
                        s1 += E1[i0 + i, k] * h1[k] * X1[i, k]; s2 += E2[i0 + i, k] * h2[k] * X2[i, k]
                    c1 = (1.0 - pi) / s1; c2 = pi / s2
                    for k in range(K):
                        P[i, k] = c1 * E1[i0 + i, k] * h1[k] * X1[i, k] + c2 * E2[i0 + i, k] * h2[k] * X2[i, k]
                for k in range(K):
                    # classify seats, then recursion over active ones
                    fixed = 0; na = 0; drop = 0.0; mu = 0.0; vs = 0.0
                    act = np.empty(nr)
                    for i in range(nr):
                        q = P[i, k]
                        mu += q; vs += q * (1.0 - q)
                        if q <= eps:
                            drop += q
                        elif q >= 1.0 - eps:
                            fixed += 1; drop += 1.0 - q
                        else:
                            act[na] = q; na += 1
                    tvr += wb * drop
                    mu_r[k] += wb * mu; m2_r[k] += wb * mu * mu
                    vseat[a, k] += wb * vs
                    d = np.zeros(na + 1); d[0] = 1.0
                    for j in range(na):
                        q = act[j]
                        for s in range(j + 1, 0, -1):
                            d[s] = d[s] * (1.0 - q) + d[s - 1] * q
                        d[0] *= 1.0 - q
                    for s in range(na + 1):
                        rd[k, s + fixed] += wb * d[s]
                    if fixed < rlo[k]:
                        rlo[k] = fixed
                    if fixed + na > rhi[k]:
                        rhi[k] = fixed + na
            tv[a] += tvr
            for k in range(K):
                mean[a, k] += mu_r[k]
                vreg[a, k] += m2_r[k] - mu_r[k] * mu_r[k]
                # convolve tot[k][lo:hi+1] with rd[k][rlo:rhi+1]
                nlo = lo[k] + rlo[k]; nhi = hi[k] + rhi[k]
                new = np.zeros(nhi - nlo + 1)
                for x in range(lo[k], hi[k] + 1):
                    tx = tot[k, x]
                    if tx == 0.0:
                        continue
                    for y in range(rlo[k], rhi[k] + 1):
                        new[x + y - nlo] += tx * rd[k, y]
                for x in range(lo[k], hi[k] + 1):
                    tot[k, x] = 0.0
                for x in range(nlo, nhi + 1):
                    tot[k, x] = new[x - nlo]
                lo[k] = nlo; hi[k] = nhi
        for k in range(K):
            for x in range(N + 1):
                dist[a, k, x] = tot[k, x]
    return dist, mean, vseat, vreg, tv


@njit(parallel=True, cache=True)
def exact_chunk_L(Dn, base, pres, Lseat, Mu, roff, rshift, rw, rcount, pi, t1, t2, eps, N):
    A, K = Dn.shape
    n = base.shape[0]
    E1 = np.zeros((n, K)); E2 = np.zeros((n, K))
    for i in range(n):
        for k in range(K):
            if pres[i, k]:
                E1[i, k] = np.exp(base[i, k] / t1); E2[i, k] = np.exp(base[i, k] / t2)
    return _exact_chunk_L(Dn, E1, E2, Lseat, Mu, roff, rshift, rw, rcount, pi, t1, t2, eps, N)
