"""v0.2 engine: outer randomised quasi-Monte Carlo over (national poll factors z, seat-type factors mu); inner:
regional Gauss-Hermite quadrature (an approximation, rule `reg_ns`) x exact conditional Poisson-binomial recursion x
convolution across regions. Seat-level sampling noise is removed by the exact conditional computation.
Error bars: the standard error estimated from independent scramblings estimates the outer RQMC sampling error under
the FIXED regional rule only; it does not include the regional-rule bias (R1 audit: Green) and, with few scramblings,
can understate rare outer contributions (R1 audit: Oth). Report regional-rule changes and a conditional RQMC
comparison alongside (scripts/reg_convergence_v02*.py, scripts/cond_mc_v02.py).
mu_df: Student-t mu draws are WINSORISED at |t| <= T_CLIP (t2: 0.985% of mass per factor beyond the cut); this is not
an unbounded t distribution."""
import time
import numpy as np
from scipy.stats import qmc, norm, chi2
from .core import exact_chunk_L
from .engine import FastEngine, KC
from .quad import tensor_gh

T_CLIP = 10.0

# Regional Gauss-Hermite rules (R1 audit, 2026-10-20). The constructor default is kept at the fast v0.2.1 rule so that
# stored v0.2/v0.2.1 runs reproduce; it has a known Green >= 4 bias of about +3e-3 against conditional RQMC references.
# Results reported in the election computation note use REG_NS_AUDITED, passed EXPLICITLY (see
# scripts/r1_note_tables_v02.py); its mean bias against regional-only RQMC is audited in
# scripts/r1_fixed_outer_reps_v02.py.
REG_NS_V021_DEFAULT = (3, 3, 2, 2, 2, 2)
REG_NS_AUDITED = (3, 3, 3, 9, 5, 3)


class QEngine2(FastEngine):
    def __init__(self, model, reg_ns=REG_NS_V021_DEFAULT, nat_on=True, mu_on=True, mu_fixed=None, eps=1e-12, chunk=256):
        super().__init__(model, nat_spec=[("gh", 1)] * model.A.shape[1], reg_ns=reg_ns, eps=eps, chunk=chunk, nat_norm=True)
        self.nat_on, self.mu_on, self.mu_fixed = nat_on, mu_on, mu_fixed
        self.L_s = model.L[:, :KC, :][self.order].copy()

    def run(self, lam=0.0, m_log2=11, reps=4, seed=0, mu_mean=None, mu_sd=None, mu_df=None, mu_scale=1.0):
        """mu_df: if set, each mu_j ~ m_j + mu_scale * s_j * t(mu_df) (independent Student-t predictive, e.g. df = n-1 = 2
        and mu_scale = sqrt(1 + 1/n) for a sample of n = 3 past transitions); implemented as a normal / chi scale mixture
        with one extra Sobol dimension per factor."""
        t0 = time.time()
        m = self.m
        J = m.L.shape[2]
        mm = m.mu_mean if mu_mean is None else np.asarray(mu_mean, float)
        ss = m.mu_sd if mu_sd is None else np.asarray(mu_sd, float)
        Dl = m.scenario_shift(lam)
        dz = m.A.shape[1] if self.nat_on else 0
        dm = J if (self.mu_on and self.mu_fixed is None) else 0
        dt = dm if (dm and mu_df is not None) else 0
        rs, rw, rc = self._reg_arrays()
        N = self.N
        dists = []; tvs = []
        for r in range(reps):
            d = dz + dm + dt
            if d > 0:
                U = qmc.Sobol(d, scramble=True, seed=seed + 7919 * r).random_base2(m_log2)
                Z = norm.ppf(np.clip(U, 1e-12, 1 - 1e-12))
            else:
                Z = np.zeros((1, 0))
            M = Z.shape[0]
            Zn = Z[:, :dz] if dz else np.zeros((M, m.A.shape[1]))
            if self.mu_fixed is not None:
                Mu = np.tile(np.asarray(self.mu_fixed, float), (M, 1))
            elif dm:
                Mu = mm + ss * Z[:, dz:dz + dm]
                if dt:
                    W = chi2.ppf(np.clip(U[:, dz + dm:dz + dm + dt], 1e-12, 1 - 1e-12), mu_df) / mu_df
                    T = np.clip(Z[:, dz:dz + dm] / np.sqrt(W), -T_CLIP, T_CLIP)   # winsorised at |t| <= 10 (t2: 0.985% per factor beyond)
                    Mu = mm + mu_scale * ss * T
            else:
                Mu = np.tile(mm, (M, 1))
            Dn = (Dl[None, :] + Zn @ m.A.T)[:, :KC].copy()
            acc = np.zeros((KC, N + 1)); tvmax = 0.0
            for s in range(0, M, self.chunk):
                sl = slice(s, s + self.chunk)
                dd, mu_, vs, vr, tv = exact_chunk_L(Dn[sl], self.base_s, self.pres_s, self.L_s, Mu[sl].copy(),
                                                     self.roff, rs, rw, rc, m.pi, m.tau1, m.tau2, self.eps, N)
                acc += dd.sum(0); tvmax = max(tvmax, float(np.max(tv)) if np.size(tv) else 0.0)
            dists.append(acc / M); tvs.append(tvmax)
        D = np.array(dists)
        # R1 review 3.2: keep the truncation audit (max over nodes of the per-node TV bound reported by the core)
        return dict(dist=D.mean(0), dist_reps=D, seconds=time.time() - t0, nodes=M * reps, lam=lam, tv_max_node=max(tvs),
                    reg_ns=self.reg_ns, eps=self.eps, m_log2=m_log2, reps=reps, seed=seed)

    def outer_points(self, r, lam=0.0, m_log2=11, seed=0, mu_mean=None, mu_sd=None, mu_df=None, mu_scale=1.0):
        """Outer points of scrambling r as used by run() with the same arguments: national shifts Dn (M, KC) and
        seat-type factor values Mu (M, J). The point construction repeats run() line by line; run() itself is left
        unchanged so that stored results reproduce. tests/test_engine2_v0_2.py checks that the two agree."""
        m = self.m
        J = m.L.shape[2]
        mm = m.mu_mean if mu_mean is None else np.asarray(mu_mean, float)
        ss = m.mu_sd if mu_sd is None else np.asarray(mu_sd, float)
        Dl = m.scenario_shift(lam)
        dz = m.A.shape[1] if self.nat_on else 0
        dm = J if (self.mu_on and self.mu_fixed is None) else 0
        dt = dm if (dm and mu_df is not None) else 0
        d = dz + dm + dt
        if d > 0:
            U = qmc.Sobol(d, scramble=True, seed=seed + 7919 * r).random_base2(m_log2)
            Z = norm.ppf(np.clip(U, 1e-12, 1 - 1e-12))
        else:
            Z = np.zeros((1, 0))
        M = Z.shape[0]
        Zn = Z[:, :dz] if dz else np.zeros((M, m.A.shape[1]))
        if self.mu_fixed is not None:
            Mu = np.tile(np.asarray(self.mu_fixed, float), (M, 1))
        elif dm:
            Mu = mm + ss * Z[:, dz:dz + dm]
            if dt:
                W = chi2.ppf(np.clip(U[:, dz + dm:dz + dm + dt], 1e-12, 1 - 1e-12), mu_df) / mu_df
                T = np.clip(Z[:, dz:dz + dm] / np.sqrt(W), -T_CLIP, T_CLIP)
                Mu = mm + mu_scale * ss * T
        else:
            Mu = np.tile(mm, (M, 1))
        Dn = (Dl[None, :] + Zn @ m.A.T)[:, :KC].copy()
        return Dn, Mu

    def seat_marginals(self, lam=0.0, m_log2=11, reps=4, seed=0, mu_mean=None, mu_sd=None, mu_df=None, mu_scale=1.0,
                       return_reps=False):
        """Seat-level win probabilities E[p_ik] under the SAME model and outer points as run() with the same arguments
        (national factors, seat-type factors mu with mu_fixed / mu_on / mu_df, regional rule reg_ns). The regional
        integral uses the Gauss-Hermite rule of run(); seat probabilities are the mixture softmax of the model.
        seat_marginals(...).sum(0) equals the mean of run(...)['dist'] up to floating point (truncation eps aside).
        Returns (N, KC) averaged over scramblings, or (reps, N, KC) if return_reps. Cost: about
        reps x 2^m_log2 x (regional nodes) x N x KC softmax evaluations in NumPy, far slower than run() for large rules.
        Replaces the inherited ExactEngine/FastEngine method, which omitted L @ mu and the outer RQMC points
        (review of 2026-10-21, issue 1)."""
        out = []
        for r in range(reps):
            Dn, Mu = self.outer_points(r, lam=lam, m_log2=m_log2, seed=seed, mu_mean=mu_mean, mu_sd=mu_sd,
                                       mu_df=mu_df, mu_scale=mu_scale)
            out.append(self._marginals_at(Dn, np.full(len(Dn), 1.0 / len(Dn)), Mu))
        out = np.array(out)
        return out if return_reps else out.mean(0)
