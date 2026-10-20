import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import sys, json, numpy as np
sys.path.insert(0, _ROOT)
from ukseat.setup2024 import make_model
from ukseat.engine import FastEngine
from ukseat.report import events
m = make_model()
NAT = [("trap", 0.5, 7.5), ("gh", 5), ("gh", 2), ("gh", 2), ("gh", 2)]
res = {}
for tag, reg in [("R0 3/2/2/2/2/2", (3, 2, 2, 2, 2, 2)), ("auto", "auto"),
                 ("auto-fine", {"hi": 0.3, "mid": 0.03, "lo": 0.003}), ("full 5/3/3/3/3/3", (5, 3, 3, 3, 3, 3))]:
    e = FastEngine(m, nat_spec=NAT, reg_ns=reg)
    if tag == "auto":
        for j, r in enumerate(m.regions):
            print(f"   {r:26s} share {np.round(e.reg_ev[j]/e.reg_ev[j].sum(),3)} -> {e.reg_ns_for(j)}")
    r = e.run(0.0)
    ev = events(r["dist"]); res[tag] = dict(ev=ev, sec=r["seconds"], nodes_reg=r["nodes_reg"], dist=r["dist"].tolist())
    print(f"{tag:18s} sec={r['seconds']:6.1f} reg nodes={r['nodes_reg']}", flush=True)
ref = res["full 5/3/3/3/3/3"]["ev"]
for tag in res:
    d = {k: res[tag]["ev"][k] - ref[k] for k in ref}
    print(f"{tag:18s} max|diff vs full| = {max(abs(v) for v in d.values()):.2e}  ({max(d, key=lambda k: abs(d[k]))})")
json.dump(res, open(f"{_ROOT}/out/reg_rule_v0_1.json", "w"))
