import numpy as np
from .model import CATS

K = {k: i for i, k in enumerate(CATS)}
# (label, party, op, threshold): ">=" -> P(S >= t), "<=" -> P(S <= t)
EVENTS = [("Lab>=326 (majority)", "Lab", ">=", 326), ("Lab>=400", "Lab", ">=", 400), ("Lab>=450", "Lab", ">=", 450),
          ("Lab<=411 (actual or fewer)", "Lab", "<=", 411), ("Lab<=250", "Lab", "<=", 250),
          ("Con>=121 (actual or more)", "Con", ">=", 121), ("Con<100", "Con", "<=", 99), ("Con<=50", "Con", "<=", 50),
          ("Con>=250", "Con", ">=", 250),
          ("LD>=72 (actual or more)", "LD", ">=", 72), ("LD>=50", "LD", ">=", 50),
          ("Ref>=5 (actual or more)", "Ref", ">=", 5), ("Ref>=50", "Ref", ">=", 50),
          ("Grn>=4 (actual or more)", "Grn", ">=", 4), ("SNP<=9 (actual or fewer)", "SNP", "<=", 9),
          ("PC>=4 (actual or more)", "PC", ">=", 4), ("Oth>=5 (actual or more)", "Oth", ">=", 5)]


def prob_from_dist(dist, party, op, t):
    d = dist[K[party]]
    return d[t:].sum() if op == ">=" else d[: t + 1].sum()


def prob_from_samples(S, party, op, t):
    x = S[:, K[party]]
    return np.mean(x >= t) if op == ">=" else np.mean(x <= t)


def summarise_dist(dist):
    out = {}
    s = np.arange(dist.shape[1])
    for k, i in K.items():
        d = dist[i]; mu = d @ s; sd = np.sqrt(d @ (s - mu) ** 2)
        c = np.cumsum(d)
        q = [int(np.searchsorted(c, x)) for x in (0.05, 0.5, 0.95)]
        out[k] = dict(mean=mu, sd=sd, q05=q[0], q50=q[1], q95=q[2])
    return out


def events(dist):
    ev = {lab: prob_from_dist(dist, p, op, t) for lab, p, op, t in EVENTS}
    ev["Hung (no party>=326)"] = 1 - sum(dist[i, 326:].sum() for i in range(dist.shape[0]))
    return ev


def events_mc(S):
    ev = {lab: prob_from_samples(S, p, op, t) for lab, p, op, t in EVENTS}
    ev["Hung (no party>=326)"] = 1 - np.mean((S >= 326).any(1))
    return ev
