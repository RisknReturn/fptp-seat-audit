"""UK Westminster seat model v0.1 (GB, 2024 boundaries).

Layers (all shocks in log-share 'delta' space):
  seat i, party k:  eta_ik = log b_ik + dbar_{B(i),k} + Delta_lambda_k + [A z]_k + [C_{r(i)} w_{r(i)}]_k
  z   ~ N(0, I_5)   national polling-error factor (Con, Lab, LD, Ref, Grn), common to GB
  w_r ~ N(0, I_dr)  regional factor, independent across regions given z
  win prob: p_ik = (1-pi) softmax_k(eta_i/tau1) + pi softmax_k(eta_i/tau2)
            (two-scale Gumbel noise on log shares; seats independent given z, w)
"""
import json
import numpy as np
import pandas as pd
from .swing import shares, solve_delta, jacobian

CATS = ["Con", "Lab", "LD", "Ref", "Grn", "SNP", "PC", "Oth"]   # contenders
ALL = CATS + ["R"]                                              # R = non-contending residual
NAT = ["Con", "Lab", "LD", "Ref", "Grn"]                        # national-factor parties
IDX = {k: i for i, k in enumerate(ALL)}
FLOOR = 0.002


def build_baseline(df):
    """2019 notional shares on 2024 boundaries with documented imputations."""
    n = len(df)
    tot = df[[c for c in df if c.startswith("n19_") and c != "n19_OthMax"]].sum(axis=1).values
    b = np.zeros((n, len(ALL)))
    b[:, IDX["Con"]] = df.n19_Con / tot
    b[:, IDX["Lab"]] = df.n19_Lab / tot
    b[:, IDX["LD"]] = np.where(df.n19_LD > 0, df.n19_LD / tot, df.s15map_LD)          # LD stood aside (2019 pact)
    b[:, IDX["Grn"]] = np.where(df.n19_Grn > 0, df.n19_Grn / tot, df.s15map_Grn)      # Greens stood aside
    b[:, IDX["Ref"]] = df.s15map_UKIP.values                                           # Reform potential = UKIP 2015 map
    b[:, IDX["SNP"]] = np.where(df.country == "Scotland", df.n19_SNP / tot, 0)
    b[:, IDX["PC"]] = np.where(df.country == "Wales", df.n19_PC / tot, 0)
    b[:, IDX["Oth"]] = df.n19_OthMax / tot
    b[:, IDX["R"]] = (df.n19_Oth - df.n19_OthMax + df.n19_UKIP) / tot
    present = np.ones_like(b, dtype=bool)
    present[:, IDX["SNP"]] = (df.country == "Scotland").values
    present[:, IDX["PC"]] = (df.country == "Wales").values
    b = np.where(present, np.maximum(b, FLOOR), 0.0)
    with np.errstate(divide="ignore"):
        logb = np.where(present, np.log(b), -np.inf)
    return logb, tot


class SeatModel:
    def __init__(self, seats_csv, calib_json, polls, sigma_nat_pts, corr_nat, reg_cfg, e_ref=None):
        df = pd.read_csv(seats_csv, index_col=0)
        self.spk_mask = (df.n19_Spk > 0).values                    # Speaker's seat held fixed
        self.n_spk = int(self.spk_mask.sum())
        df = df[~self.spk_mask].copy()
        self.df = df
        self.n = len(df)
        self.logb, self.w = build_baseline(df)
        wcols = [f"n19_{k}" for k in ["Con", "Lab", "LD", "Ref", "Grn", "SNP", "PC", "OthMax"]]
        self.n19_winner = df[wcols].astype(float).idxmax(axis=1).str[4:].replace({"OthMax": "Oth"}).values
        if getattr(self, "_cov_csv", None):
            from .model2 import build_baseline2
            cov = pd.read_csv(self._cov_csv, index_col=0).loc[df.index]
            self.logb = build_baseline2(self, cov)
        self.block = df.country.map({"England": 0, "Scotland": 1, "Wales": 2}).values
        self.regions = sorted(df.region.unique())
        self.reg = df.region.map({r: j for j, r in enumerate(self.regions)}).values
        cal = json.load(open(calib_json))
        self.pi, self.tau1, self.tau2 = cal["mix_pi"], cal["mix_tau1"], cal["mix_tau2"]
        self.polls = polls
        self.e_ref = e_ref
        # ---- block means (raking to block targets)
        self.dbar = np.zeros((3, len(ALL)))
        tgts = self._block_targets(polls)
        for B in range(3):
            m = self.block == B
            lb = np.where(np.isfinite(self.logb[m]), self.logb[m], -60.0)
            free = np.zeros(len(ALL), bool)
            for k in tgts[B]:
                free[IDX[k]] = True
            t = np.zeros(len(ALL))
            for k, v in tgts[B].items():
                t[IDX[k]] = v
            self.dbar[B] = solve_delta(lb, self.w[m], t, free)
        # ---- national factor loading A (delta-space) from points covariance via GB Jacobian
        sd = np.array([sigma_nat_pts[k] for k in NAT]) / 100.0
        Sig = np.outer(sd, sd) * corr_nat
        J = self.gb_jacobian(np.zeros(len(ALL)))[np.ix_([IDX[k] for k in NAT], [IDX[k] for k in NAT])]
        self.Sig_nat_pts = Sig
        L = np.linalg.cholesky(Sig)
        A5 = np.linalg.solve(J, L)                                   # 5x5: delta_NAT = A5 z
        self.A = np.zeros((len(ALL), len(NAT)))
        self.A[[IDX[k] for k in NAT], :] = A5
        # ---- regional factors: covariance in delta space per region
        self.reg_cfg = reg_cfg
        self.C = []
        for j, r in enumerate(self.regions):
            parties = reg_cfg["parties_scot"] if r == "Scotland" else reg_cfg["parties_wales"] if r == "Wales" else reg_cfg["parties_eng"]
            S = np.array(reg_cfg["cov"][r if r in reg_cfg["cov"] else "England"])
            C = np.zeros((len(ALL), len(parties)))
            # eigen-factorisation so callers can truncate by variance if they wish
            ev, U = np.linalg.eigh(S)
            order = np.argsort(ev)[::-1]
            ev, U = ev[order], U[:, order]
            for a, k in enumerate(parties):
                C[IDX[k], :] = U[a, :] * np.sqrt(np.maximum(ev, 0))
            self.C.append(C)
        self.Delta = {}

    # ------------------------------------------------------------------ helpers
    def _block_targets(self, polls):
        """polls: dict with 'GB', 'Scotland', 'Wales' share dicts (percent). England is implied."""
        wS = self.w[self.block == 1].sum() / self.w.sum()
        wW = self.w[self.block == 2].sum() / self.w.sum()
        wE = 1 - wS - wW
        S, W, G = polls["Scotland"], polls["Wales"], polls["GB"]
        E = {k: (G[k] - wS * S.get(k, 0) - wW * W.get(k, 0)) / wE for k in NAT}
        return [{k: E[k] / 100 for k in NAT},
                {k: S[k] / 100 for k in NAT + ["SNP"]},
                {k: W[k] / 100 for k in NAT + ["PC"]}]

    def mean_delta(self):
        """(n, K) seat-level delta at the national/regional centre (no shocks)."""
        return self.dbar[self.block]

    def gb_shares(self, Delta):
        d = self.mean_delta() + Delta
        lb = np.where(np.isfinite(self.logb), self.logb, -60.0)
        v = shares(lb, d)
        return (self.w[:, None] * v).sum(0) / self.w.sum()

    def gb_jacobian(self, Delta, h=1e-6):
        K = len(ALL)
        J = np.zeros((K, K))
        f0 = self.gb_shares(Delta)
        for j in range(K):
            e = np.zeros(K); e[j] = h
            J[:, j] = (self.gb_shares(Delta + e) - f0) / h
        return J

    def scenario_shift(self, lam):
        """Common delta shift that moves GB shares by -lam * e_ref (poll minus result, points)."""
        if lam == 0 or self.e_ref is None:
            return np.zeros(len(ALL))
        if lam in self.Delta:
            return self.Delta[lam]
        target = self.gb_shares(np.zeros(len(ALL))).copy()
        for k in NAT:
            target[IDX[k]] -= lam * self.e_ref[k] / 100.0
        D = np.zeros(len(ALL))
        fi = [IDX[k] for k in NAT]
        for _ in range(50):
            g = self.gb_shares(D)[fi] - target[fi]
            if np.max(np.abs(g)) < 1e-12:
                break
            J = self.gb_jacobian(D)[np.ix_(fi, fi)]
            D[fi] -= np.linalg.solve(J, g)
        self.Delta[lam] = D
        return D

    def win_probs(self, eta, present):
        """eta: (..., n, K) with -inf for absent parties. Returns (..., n, K)."""
        out = 0.0
        for wgt, tau in ((1 - self.pi, self.tau1), (self.pi, self.tau2)):
            x = np.where(present, eta / tau, -np.inf)
            m = np.max(x, axis=-1, keepdims=True)
            e = np.exp(x - m)
            out = out + wgt * e / e.sum(axis=-1, keepdims=True)
        return out
