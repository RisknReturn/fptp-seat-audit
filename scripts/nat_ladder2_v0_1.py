"""One-at-a-time refinement of every national dimension (raw vs party-normalised rotation), regional R1."""
import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import sys, json, numpy as np
sys.path.insert(0, _ROOT)
from ukseat.setup2024 import make_model
from ukseat.engine import FastEngine
from ukseat.report import events
m = make_model()
mc = json.load(open(f"{_ROOT}/out/mc_compare_std_v0_1.json"))
m6 = [r for r in mc["runs"] if r["N"] == 1_000_000]
mcev = {k: float(np.mean([r["ev"][k] for r in m6])) for k in m6[0]["ev"]}
R1 = (3, 3, 2, 2, 2, 2)
base = [("trap", 0.5, 7.5), ("gh", 5), ("gh", 3), ("gh", 5), ("gh", 3)]
ref = {0: ("trap", 0.35, 7.5), 1: ("gh", 7), 2: ("gh", 5), 3: ("gh", 7), 4: ("gh", 5)}
rows = []
for norm in (False, True):
    e0 = FastEngine(m, nat_spec=base, reg_ns=R1, nat_norm=norm)
    print("norm", norm, "impact shares", np.round(e0.nat_ev / e0.nat_ev.sum(), 3), flush=True)
    for tag, spec in [("base", base)] + [(f"refine y{j+1}", base[:j] + [ref[j]] + base[j + 1:]) for j in range(5)]:
        e = FastEngine(m, nat_spec=spec, reg_ns=R1, nat_norm=norm); r = e.run(0.0); ev = events(r["dist"])
        rows.append(dict(norm=norm, tag=tag, nodes=r["nodes_outer"], sec=r["seconds"], ev=ev))
        dmc = {k: ev[k] - mcev[k] for k in ev}; km = max(dmc, key=lambda k: abs(dmc[k]))
        print(f"  {tag:10s} nodes={r['nodes_outer']:6d} sec={r['seconds']:5.1f} max|ev-MC|={abs(dmc[km]):.4f} ({km})", flush=True)
json.dump(dict(rows=rows, mc=mcev), open(f"{_ROOT}/out/nat_ladder2_v0_1.json", "w"))
