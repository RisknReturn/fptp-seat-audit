"""Reference exact run (lambda=0): N5 national x R3 regional."""
import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import sys, json
sys.path.insert(0, _ROOT)
from ukseat.setup2024 import make_model
from ukseat.engine import FastEngine
from ukseat.report import events, summarise_dist
m = make_model()
NAT = [("trap", 0.25, 7.5), ("gh", 9), ("gh", 5), ("gh", 3), ("gh", 3)]
r = FastEngine(m, nat_spec=NAT, reg_ns=(5, 3, 3, 3, 3, 3)).run(0.0)
json.dump(dict(ev=events(r["dist"]), summ=summarise_dist(r["dist"]), sec=r["seconds"], nodes=r["nodes_outer"],
               reg_nodes=r["nodes_reg"], tv=r["tv_bound"], dist=r["dist"].tolist(),
               var=dict(nat=r["var_nat"].tolist(), reg=r["var_reg"].tolist(), seat=r["var_seat"].tolist())),
          open(f"{_ROOT}/out/ref_exact_v0_1.json", "w"))
print("ref exact sec %.0f" % r["seconds"])
