"""Proportional (log-ratio) swing with raking to block-level targets.

shares v_ik(delta) = b_ik * exp(delta_k) / sum_j b_ij * exp(delta_j)
delta is solved so that the weight-averaged shares match target totals.
"""
import numpy as np


def shares(logb, delta):
    """logb: (n, K) log baseline shares (-inf = party absent); delta: (K,) or (..., K)."""
    x = logb + delta
    m = np.max(x, axis=-1, keepdims=True)
    e = np.exp(x - m)
    return e / e.sum(axis=-1, keepdims=True)


def block_totals(logb, delta, w):
    v = shares(logb, delta)
    return (w[:, None] * v).sum(0) / w.sum()


def solve_delta(logb, w, target, free, delta0=None, tol=1e-12, it=100):
    """Find delta (with delta[~free] fixed) so block totals on `free` parties hit target[free].
    Newton on the free coordinates. target: (K,) shares (only free entries used)."""
    K = logb.shape[1]
    d = np.zeros(K) if delta0 is None else delta0.copy()
    fi = np.where(free)[0]
    for _ in range(it):
        v = shares(logb, d)
        W = w / w.sum()
        tot = (W[:, None] * v).sum(0)
        g = tot[fi] - target[fi]
        if np.max(np.abs(g)) < tol:
            return d
        # d tot_k / d delta_j = sum_i W_i v_ik (1[k=j] - v_ij)
        Jm = -(W[:, None, None] * v[:, :, None] * v[:, None, :]).sum(0)
        Jm[np.diag_indices(K)] += tot
        step = np.linalg.lstsq(Jm[np.ix_(fi, fi)], g, rcond=None)[0]
        # damped Newton: backtrack until the gap shrinks
        f0 = np.sum(g ** 2); a = 1.0
        while a > 1e-6:
            dn = d.copy(); dn[fi] -= a * np.clip(step, -5, 5)
            tn = (W[:, None] * shares(logb, dn)).sum(0)
            if np.sum((tn[fi] - target[fi]) ** 2) < f0:
                break
            a *= 0.5
        d = dn
    raise RuntimeError("raking did not converge: max gap %.2e" % np.max(np.abs(g)))


def jacobian(logb, w, delta):
    """d(block totals)/d(delta) at delta: (K, K)."""
    W = w / w.sum()
    v = shares(logb, delta)
    tot = (W[:, None] * v).sum(0)
    Jm = -(W[:, None, None] * v[:, :, None] * v[:, None, :]).sum(0)
    Jm[np.diag_indices(len(tot))] += tot
    return Jm
