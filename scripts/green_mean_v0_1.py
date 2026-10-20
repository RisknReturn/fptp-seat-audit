"""Where the remaining Green error comes from: E[Green seats] = sum_i E[p_i,Grn] (no recursion needed),
under different allocations of quadrature points. Target: plain MC mean over 3e6 draws."""
import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import sys, json, time, numpy as np
sys.path.insert(0, _ROOT)
from ukseat.setup2024 import make_model
from ukseat.engine import ExactEngine
from ukseat.model import CATS
m = make_model()
mc = json.load(open(f"{_ROOT}/out/mc_compare_std_v0_1.json"))
m6 = [r for r in mc["runs"] if r["N"] == 1_000_000]
g = CATS.index("Grn")
mcm = np.mean([r["mean"][g] for r in m6]); mcse = np.mean([r["sd"][g] for r in m6]) / np.sqrt(3e6)
print(f"MC 3e6 E[Grn] = {mcm:.4f} (se {mcse:.4f})")
class E2(ExactEngine):
    nat_norm = True
rows = []
for nat in [(3, 3, 3, 3, 3), (3, 3, 3, 3, 7)]:
    for reg in [(5, 5, 3, 2, 2, 2), (5, 3, 3, 3, 3, 3), (3, 3, 3, 3, 5, 5), (3, 3, 3, 5, 5, 5)]:
        t = time.time(); e = E2(m, nat_ns=nat, reg_ns=reg)
        P = e.seat_marginals(0.0); v = P[:, g].sum()
        rows.append(dict(nat=nat, reg=reg, mean=float(v), z=float((v - mcm) / mcse), sec=time.time() - t))
        print(f"nat {nat} reg {reg}: E[Grn]={v:.4f}  z={((v-mcm)/mcse):+.1f}  ({time.time()-t:.0f}s)", flush=True)
json.dump(dict(mc=mcm, se=mcse, rows=rows), open(f"{_ROOT}/out/green_mean_v0_1.json", "w"))
