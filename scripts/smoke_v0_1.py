import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import sys, time, numpy as np
sys.path.insert(0, _ROOT)
from ukseat.setup2024 import make_model
from ukseat.engine import ExactEngine, monte_carlo
from ukseat.model import CATS
m = make_model()
print("A rows (Con..Grn):\n", np.round(m.A[:5], 3))
for ns, rs in [((3, 3, 2, 2, 2), (3, 2, 2, 2, 2, 2)), ((5, 3, 3, 2, 2), (3, 3, 2, 2, 2, 2))]:
    e = ExactEngine(m, nat_ns=ns, reg_ns=rs)
    r = e.run(0.0)
    print(ns, rs, "outer nodes", r["nodes_outer"], "sec %.1f" % r["seconds"], "tv %.1e" % r["tv_bound"])
    print("  mean", dict(zip(CATS, np.round(r["mean"], 1))))
    sd = np.sqrt(r["var_seat"] + r["var_reg"] + r["var_nat"])
    print("  sd  ", dict(zip(CATS, np.round(sd, 1))))
    print("  var shares nat/reg/seat (Lab):", np.round(np.array([r["var_nat"][1], r["var_reg"][1], r["var_seat"][1]]) / sd[1] ** 2, 3))
    print("  P(Lab>=326) %.6f  P(Lab<=411) %.6f  sum dist %.12f" % (r["dist"][1, 326:].sum(), r["dist"][1, :412].sum(), r["dist"][1].sum()))
out, t = monte_carlo(m, 0.0, 20000)
print("MC 20k sec %.1f mean" % t, dict(zip(CATS, out.mean(0).round(1))), "sd", dict(zip(CATS, out.std(0).round(1))))
print("MC P(Lab<=411) %.4f" % np.mean(out[:, 1] <= 411))
