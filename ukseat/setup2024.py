import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import numpy as np
from .model import SeatModel
from . import config_2024 as C

ROOT = _ROOT
ENG = ["Con", "Lab", "LD", "Ref", "Grn"]


def reg_cfg():
    import json
    cal = json.load(open(f"{ROOT}/data/calib_v0_1.json"))
    S3 = np.array(cal["cov_reg_ConLabLD"])
    S = np.zeros((5, 5)); S[:3, :3] = S3; S[3, 3] = 0.20 ** 2; S[4, 4] = 0.25 ** 2
    Ss = np.zeros((6, 6)); Ss[:5, :5] = S; Ss[5, 5] = 0.08 ** 2     # SNP (Scottish poll error)
    Sw = np.zeros((6, 6)); Sw[:5, :5] = S; Sw[5, 5] = 0.12 ** 2     # Plaid (Welsh poll error)
    return {"parties_eng": ENG, "parties_scot": ENG + ["SNP"], "parties_wales": ENG + ["PC"],
            "cov": {"England": S.tolist(), "Scotland": Ss.tolist(), "Wales": Sw.tolist()}}


def make_model(polls=None, e_ref=None, sigma=None):
    return SeatModel(f"{ROOT}/data/gb2024_seats_v0_1.csv", f"{ROOT}/data/calib_v0_1.json",
                     polls or C.POLLS_2024, sigma or C.SIGMA_NAT, C.CORR_NAT, reg_cfg(),
                     e_ref=C.E_REF_2015 if e_ref is None else e_ref)
