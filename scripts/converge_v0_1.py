"""Quadrature ladders (national rule at fixed regional rule, and regional rule at fixed national rule)."""
import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import sys, json, numpy as np
sys.path.insert(0, _ROOT)
from ukseat.setup2024 import make_model
from ukseat.engine import FastEngine
from ukseat.report import events, summarise_dist

m = make_model()
OUT = f"{_ROOT}/out/converge_v0_1.json"
R0 = (3, 2, 2, 2, 2, 2)
NATS = {"N1": [("gh", 3), ("gh", 2), ("gh", 2), ("gh", 2), ("gh", 2)],
        "N2": [("trap", 0.5, 7.5), ("gh", 3), ("gh", 2), ("gh", 2), ("gh", 2)],
        "N3": [("trap", 0.35, 7.5), ("gh", 5), ("gh", 2), ("gh", 2), ("gh", 2)],
        "N4": [("trap", 0.35, 7.5), ("gh", 7), ("gh", 3), ("gh", 2), ("gh", 2)],
        "N5": [("trap", 0.25, 7.5), ("gh", 9), ("gh", 5), ("gh", 3), ("gh", 3)],
        "N6": [("trap", 0.25, 7.5), ("gh", 11), ("gh", 7), ("gh", 3), ("gh", 3)]}
REGS = {"R0": R0, "R1": (3, 3, 2, 2, 2, 2), "R2": (5, 3, 3, 2, 2, 2), "R3": (5, 3, 3, 3, 3, 3),
        "R4": (5, 5, 3, 3, 3, 3), "R5": (7, 5, 3, 3, 3, 3)}
rows = []
FastEngine(m, nat_spec=NATS["N1"], reg_ns=R0).run(0.0)  # compile
def run(ladder, nk, rk):
    e = FastEngine(m, nat_spec=NATS[nk], reg_ns=REGS[rk])
    r = e.run(0.0)
    rows.append(dict(ladder=ladder, nat=nk, reg=rk, nodes=r["nodes_outer"], reg_nodes=r["nodes_reg"], sec=r["seconds"],
                     ev=events(r["dist"]), summ=summarise_dist(r["dist"]), tv=r["tv_bound"]))
    json.dump(rows, open(OUT, "w"))
    print(f"{ladder} {nk} {rk} nodes={r['nodes_outer']} sec={r['seconds']:.1f}", flush=True)
for nk in NATS:
    run("nat", nk, "R0")
for rk in REGS:
    run("reg", "N2", rk)
