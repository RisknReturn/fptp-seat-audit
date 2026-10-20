"""v0.2 seat model: v0.1 + (i) 2024 nominations, (ii) seat-type factors mu with seat-specific loadings.

eta_ik = log b_ik + dbar_Bk + Delta_k + [A z]_k + [C_r w_r]_k + sum_j L_ikj mu_j
  mu_T  anti-Conservative tactical consolidation: Con-held (2019 notional) seats, +1 on the larger of Lab/LD, -1 on the other
  mu_C  beyond-proportional Conservative change: (b_Con - regional mean b_Con) on Con
  mu_R  Leave -> Reform: (Leave - 52)/10 on Ref
  mu_L  Muslim share -> Labour: (Muslim - mean)/10 on Lab, and -KAPPA*(Muslim/10) on 'Oth' where an independent or
        Workers Party candidate stands (one factor moves both; KAPPA is a scenario parameter, see doc)
Loadings are centred per party (vote-weighted) so a factor redistributes votes across seats without moving the
national share to first order. mu ~ N(m, s) with m, s from 2010-2019 transitions unless a scenario fixes them.
"""
import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import json
import numpy as np
import pandas as pd
from .model import SeatModel, ALL, CATS, IDX, FLOOR

FACTORS = ["mu_T", "mu_C", "mu_R", "mu_L"]
IND_MP_RATIO = 0.33   # sitting MP standing as independent: share = 0.33 x 2019 party share
                      # (2019 precedents: Field 0.22, Grieve 0.50, Gauke 0.46, Soubry 0.15 of their 2017 party share)
OTH_DEFAULT = 0.02    # baseline for an independent / minor-party contender with no history in the seat


class SeatModel2(SeatModel):
    def __init__(self, *a, cov_csv=None, kappa=4.0, mu_mean=None, mu_sd=None, neutral=True, **k):
        self._cov_csv, self.kappa, self.neutral = cov_csv, kappa, neutral
        super().__init__(*a, **k)
        cov = pd.read_csv(cov_csv, index_col=0).loc[self.df.index]
        self.cov = cov
        n = self.n
        hist = json.load(open(f"{_ROOT}/data/calib_mu_v02.json"))["summary"]
        self.mu_mean = np.array([hist[f]["mean"] for f in FACTORS]) if mu_mean is None else np.asarray(mu_mean, float)
        self.mu_sd = np.array([hist[f]["sd"] for f in FACTORS]) if mu_sd is None else np.asarray(mu_sd, float)
        # loadings (n, K_all, J)
        L = np.zeros((n, len(ALL), len(FACTORS)))
        b = np.exp(np.where(np.isfinite(self.logb), self.logb, -60))
        con_held = self.n19_winner == "Con"
        lab_ch = b[:, IDX["Lab"]] >= b[:, IDX["LD"]]
        ew = (self.df.country != "Scotland").values
        sel = con_held & ew
        L[sel & lab_ch, IDX["Lab"], 0] = 1; L[sel & lab_ch, IDX["LD"], 0] = -1
        L[sel & ~lab_ch, IDX["LD"], 0] = 1; L[sel & ~lab_ch, IDX["Lab"], 0] = -1
        bc = b[:, IDX["Con"]]
        L[:, IDX["Con"], 1] = bc - pd.Series(bc, index=self.df.index).groupby(self.df.region.values).transform("mean").values
        L[:, IDX["Ref"], 2] = (cov.leave.values - 52) / 10
        mus = cov.muslim.values
        L[:, IDX["Lab"], 3] = (mus - mus.mean()) / 10
        has_minor = (cov.wpb24.values | (cov.stand24_Oth.values & ~cov.ind_sitting_mp.values))
        L[has_minor, IDX["Oth"], 3] = -kappa * mus[has_minor] / 10
        # centre per party with vote weights at the centre shares
        v = self.centre_shares()
        for j in range(len(FACTORS)):
            for kk in range(len(ALL)):
                wv = self.w * v[:, kk]
                if wv.sum() > 0 and np.any(L[:, kk, j] != 0):
                    L[:, kk, j] -= (wv @ L[:, kk, j]) / wv.sum()
        L[~np.isfinite(self.logb)] = 0.0
        self.L = L
        if neutral:
            self._neutralise()

    def rerake(self, mu):
        """Exact block re-raking at a fixed mu (peer review 1.2): block deltas such that every targeted block share
        holds with L.mu included. Returns the new (3, K) dbar; the caller sets self.dbar (and restores it)."""
        from .swing import solve_delta
        tgts = self._block_targets(self.polls)
        new = self.dbar.copy()
        for B in range(3):
            m = self.block == B
            lb = np.where(np.isfinite(self.logb[m]), self.logb[m], -60.0) + self.L[m] @ np.asarray(mu, float)
            free = np.zeros(len(ALL), bool); t = np.zeros(len(ALL))
            for k, v in tgts[B].items():
                free[IDX[k]] = True; t[IDX[k]] = v
            new[B] = solve_delta(lb, self.w[m], t, free, delta0=self.dbar[B])
        return new

    def _neutralise(self, h=1e-5):
        """v0.2.1 (peer review 1.2): per-party centring does not keep vote shares fixed after the softmax, because
        dV_k/dmu = sum_i w_i p_ik (L_ik - sum_l p_il L_il). Add block x party constants c_Bk(j) to the loadings so that
        dV_B/dmu_j = 0 at mu = 0 in every block (England / Scotland / Wales): solve J_B c = -g_B with J_B the block
        share Jacobian w.r.t. a common delta. Exact to first order; second-order drift is reported by the check script."""
        from .swing import shares
        lb = np.where(np.isfinite(self.logb), self.logb, -60.0); d0 = self.mean_delta(); K = len(ALL)
        for B in np.unique(self.block):
            m = self.block == B; w = self.w[m] / self.w[m].sum(); lbB, dB = lb[m], d0[m]
            V = lambda dd: w @ shares(lbB, dB + dd)
            J = np.column_stack([(V(np.eye(K)[k] * h) - V(-np.eye(K)[k] * h)) / (2 * h) for k in range(K)])
            v0 = V(np.zeros(K))
            ref = int(np.argmax(v0))                                   # shares sum to one: fix the largest party
            fr = np.where((v0 > 0.005) & (np.arange(K) != ref))[0]      # parties with >= 0.5% in the block
            for j in range(self.L.shape[2]):
                LB = self.L[m, :, j]
                g = (w @ shares(lbB, dB + h * LB) - w @ shares(lbB, dB - h * LB)) / (2 * h)
                c = np.zeros(K)
                c[fr] = np.linalg.solve(J[np.ix_(fr, fr)], -g[fr])       # unique; no common-mode or tiny-party blow-up
                self.L[m, :, j] += c[None, :]
        self.L[~np.isfinite(self.logb)] = 0.0

    def centre_shares(self):
        from .swing import shares
        lb = np.where(np.isfinite(self.logb), self.logb, -60.0)
        return shares(lb, self.mean_delta())


def build_baseline2(model, cov):
    """Modify the v0.1 baseline with 2024 nominations (pre-election information)."""
    df = model.df
    logb = model.logb.copy()
    # parties not nominated are absent
    for k in ["Con", "Lab", "LD", "Ref", "Grn", "SNP", "PC"]:
        absent = ~cov[f"stand24_{k}"].values.astype(bool)
        logb[absent, IDX[k]] = -np.inf
    # 'Oth' contender
    oth = np.full(len(df), -np.inf)
    ret = cov.oth_returning_share.values
    oth[ret > 0] = np.log(ret[ret > 0])
    tot = df[[c for c in df if c.startswith("n19_") and c != "n19_OthMax"]].sum(axis=1).values
    for i in np.where(cov.ind_sitting_mp.values)[0]:
        # former party = 2019 notional winner (the sitting MP's party)
        winner = df.iloc[i][[f"n19_{k}" for k in ["Con", "Lab", "LD", "SNP", "PC", "Grn"]]].astype(float).idxmax()
        share = df.iloc[i][winner] / tot[i]
        oth[i] = np.log(max(0.33 * share, FLOOR))
    minor = (cov.wpb24.values | cov.stand24_Oth.values.astype(bool)) & ~np.isfinite(oth)
    oth[minor] = np.log(OTH_DEFAULT)
    logb[:, IDX["Oth"]] = oth
    return logb
