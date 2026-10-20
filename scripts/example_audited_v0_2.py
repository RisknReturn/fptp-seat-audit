"""Minimal example (v0.2): event probabilities of the fixed v0.2.1 model with the AUDITED regional rule.

    python scripts/example_audited_v0_2.py [m_log2] [reps]

Defaults m_log2 = 8, reps = 4 (1,024 outer points; a few minutes on two cores). The note's tables use 2^11 x 4 and
larger designs. The engine's constructor default regional rule (3,3,2,2,2,2) is kept only so that stored v0.2/v0.2.1
runs reproduce; it has a known Green >= 4 bias of about +3e-3, so pass REG_NS_AUDITED explicitly as below."""
import os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import numpy as np
from ukseat.model2 import SeatModel2
from ukseat.setup2024 import reg_cfg
from ukseat import config_2024 as C
from ukseat.engine2 import QEngine2, REG_NS_AUDITED
from ukseat.report import events

m_log2 = int(sys.argv[1]) if len(sys.argv) > 1 else 8
reps = int(sys.argv[2]) if len(sys.argv) > 2 else 4
model = SeatModel2(f"{ROOT}/data/gb2024_seats_v0_1.csv", f"{ROOT}/data/calib_v0_1.json", C.POLLS_2024, C.SIGMA_NAT,
                   C.CORR_NAT, reg_cfg(), e_ref=C.E_REF_2015, cov_csv=f"{ROOT}/data/gb2024_v02_covariates.csv",
                   kappa=4.0, neutral=True)
r = QEngine2(model, reg_ns=REG_NS_AUDITED).run(lam=0.0, m_log2=m_log2, reps=reps, seed=0)
E = np.array([[v for v in events(d).values()] for d in r["dist_reps"]])
names = list(events(r["dist"]).keys())
print(f"regional rule {r['reg_ns']}, outer 2^{m_log2} x {reps}, {r['seconds']:.0f} s")
for i, k in enumerate(names):
    label = "No category >= 326 (Oth counted as one category)" if k.startswith("Hung") else k
    print(f"  {label:48s} {E[:, i].mean():.5f}  (scramble SE {E[:, i].std(ddof=1) / np.sqrt(reps):.1e})")
