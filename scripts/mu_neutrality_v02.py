"""Reviewer check 1.2: does per-party centring of mu loadings keep national (GB) and block vote shares fixed
after the softmax?  Finite differences at mu=0 and the realised 2024 mu, oracle model (actual shares)."""
import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import sys, json, numpy as np
sys.path.insert(0, _ROOT)
sys.argv = [sys.argv[0], "4"]
exec(open(f"{_ROOT}/scripts/v02_experiments.py").read().split("out = {")[0])
from ukseat.swing import shares
from ukseat.model import ALL
mo = mk(POLLS_ACT)
lb = np.where(np.isfinite(mo.logb), mo.logb, -60.0); d0 = mo.mean_delta()
def nat(mu, mask=None):
    v = shares(lb, d0 + np.einsum("ikj,j->ik", mo.L, mu))
    w = mo.w if mask is None else mo.w * mask
    return (w[:, None] * v).sum(0) / w.sum() * 100
show = ["Con", "Lab", "LD", "Ref", "Grn", "SNP", "PC", "Oth"]
base = nat(np.zeros(4)); out = {"base_GB": dict(zip(ALL, base.round(3)))}
print("GB shares at mu=0:", {k: round(base[ALL.index(k)], 2) for k in show})
for j, f in enumerate(FACTORS):
    e = np.zeros(4); e[j] = 1e-4
    g = (nat(e) - nat(-e)) / 2e-4
    mf = np.zeros(4); mf[j] = MU_REAL[j]
    dr = nat(mf) - base
    out[f] = dict(dV_dmu=dict(zip(ALL, g.round(3))), shift_at_realised=dict(zip(ALL, dr.round(3))))
    print(f"{f}: dV/dmu (pt per unit) " + " ".join(f"{k} {g[ALL.index(k)]:+.2f}" for k in show))
    print(f"      shift at realised {MU_REAL[j]:+.3f}: " + " ".join(f"{k} {dr[ALL.index(k)]:+.2f}" for k in show))
dr = nat(MU_REAL) - base
print("all four realised: " + " ".join(f"{k} {dr[ALL.index(k)]:+.2f}" for k in show))
for B, nm in enumerate(["England", "Scotland", "Wales"]):
    msk = (mo.block == B).astype(float)
    db = nat(MU_REAL, msk) - nat(np.zeros(4), msk)
    print(f"   {nm}: " + " ".join(f"{k} {db[ALL.index(k)]:+.2f}" for k in show))
out["all_realised_shift"] = dict(zip(ALL, dr.round(3)))
json.dump(out, open(f"{ROOT}/out/mu_neutrality_v02.json", "w"), indent=1)
