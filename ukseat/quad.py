"""Gaussian quadrature rules for N(0, I_d): anisotropic tensor Gauss-Hermite."""
import itertools
import numpy as np


def gh1(n):
    """Probabilists' Gauss-Hermite nodes/weights for N(0,1)."""
    if n == 1:
        return np.zeros(1), np.ones(1)
    x, w = np.polynomial.hermite_e.hermegauss(n)
    return x, w / w.sum()


def tensor_gh(ns):
    """ns: points per dimension. Returns nodes (M, d), weights (M,)."""
    rules = [gh1(n) for n in ns]
    X = np.array(list(itertools.product(*[r[0] for r in rules])))
    W = np.array([np.prod(c) for c in itertools.product(*[r[1] for r in rules])])
    return X.reshape(len(W), len(ns)), W
