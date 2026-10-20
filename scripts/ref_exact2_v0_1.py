"""Reference exact run v2 (lambda=0): party-normalised national rotation, every dimension refined."""
import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import sys, json
sys.path.insert(0, _ROOT)
from ukseat.setup2024 import make_model
from ukseat.engine import FastEngine
from ukseat.report import events, summarise_dist
m = make_model()
NAT = [("trap", 0.35, 7.5), ("gh", 7), ("gh", 5), ("gh", 5), ("gh", 3)]
REG = (5, 5, 3, 2, 2, 2)
r = FastEngine(m, nat_spec=NAT, reg_ns=REG, nat_norm=True).run(0.0)
json.dump(dict(ev=events(r["dist"]), summ=summarise_dist(r["dist"]), sec=r["seconds"], nodes=r["nodes_outer"],
               reg_nodes=r["nodes_reg"], tv=r["tv_bound"], dist=r["dist"].tolist(), nat=str(NAT), reg=str(REG),
               var=dict(nat=r["var_nat"].tolist(), reg=r["var_reg"].tolist(), seat=r["var_seat"].tolist())),
          open(f"{_ROOT}/out/ref_exact2_v0_1.json", "w"))
print("ref exact2 sec %.0f nodes %d reg %s" % (r["seconds"], r["nodes_outer"], r["nodes_reg"]))
