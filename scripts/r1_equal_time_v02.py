"""R1 review 4.2: equal-time accuracy of the hybrid engine vs conditional RQMC / conditional MC.

Reference: large conditional RQMC run (no regional quadrature, unbiased up to its scramble SE).
Per-lambda cost: the hybrid computes one lambda per run; cond_mc computes lambda = 0 and 0.5 in the same pass, so
its per-lambda cost is taken as half of its measured time (the Poisson-binomial recursion dominates both passes).
Timings were measured on 2 shared cores and are indicative only.
For a budget B (the cost of one default hybrid run, 2^11 x 4):
  hybrid  : error = (hybrid mean - reference) [rule bias, measured] with scramble SE from its 4 replicates
  cond RQMC / cond iid : unbiased; SE = per-replicate sd / sqrt(B / cost per replicate)
RMSE_hybrid = sqrt(bias^2 + SE^2), with the reference SE reported alongside (bias is only resolved where |bias| is
well above it). Writes out/r1_equal_time_v02.json.
"""
import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import json, numpy as np

ref = json.load(open(f"{_ROOT}/out/cond_rqmc_N8192_R32_s20261022.json"))
cq = json.load(open(f"{_ROOT}/out/cond_rqmc_N2048_R16_s20261020.json"))
ci = json.load(open(f"{_ROOT}/out/cond_iid_N2048_R16_s20261021.json"))
rg = json.load(open(f"{_ROOT}/out/reg_convergence_v02b.json"))
hl = json.load(open(f"{_ROOT}/out/hybrid_lam_crn_v02.json"))
EVK = ref["events"]
R = np.array(ref["reps_events"]["0.0"]); ref_m, ref_se = R.mean(0), R.std(0, ddof=1) / np.sqrt(len(R))
B = float(np.mean([rg["rules"]["(3, 3, 2, 2, 2, 2)"]["sec"], rg["paired"]["(3, 3, 2, 2, 2, 2)"]["sec"]]))   # uncontended runs
hyb = {}
for rule in ("(3, 3, 2, 2, 2, 2)", "(3, 3, 3, 3, 3, 3)", "(3, 3, 3, 5, 5, 5)"):
    E = np.array(rg["paired"][rule]["reps"])
    hyb[rule] = dict(mean=E.mean(0), se=E.std(0, ddof=1) / 2, sec=rg["paired"][rule]["sec"])
out = dict(budget_s=B, reference=dict(N=ref["N"], reps=ref["reps"], seconds=ref["total_seconds"]), rows=[])
for name, d in (("cond_rqmc", cq), ("cond_iid", ci)):
    out[f"{name}_sec_per_rep_per_lambda"] = float(np.mean(d["seconds"]) / 2)
for k, e in enumerate(EVK):
    row = dict(event=e, ref=float(ref_m[k]), ref_se=float(ref_se[k]))
    for rule, h in hyb.items():
        bias = float(h["mean"][k] - ref_m[k])
        scale = B / h["sec"]                       # replicates of this rule affordable in B (relative to 4)
        se_b = float(h["se"][k] / np.sqrt(max(scale, 1e-9)))
        row[f"hyb{rule.replace(' ', '')}"] = dict(bias=bias, bias_z=bias / float(np.hypot(ref_se[k], h["se"][k])),
                                                   se_at_B=se_b, rmse_at_B=float(np.hypot(bias, se_b)), sec=h["sec"])
    for name, d in (("cond_rqmc", cq), ("cond_iid", ci)):
        E = np.array(d["reps_events"]["0.0"])
        n_rep = B / (np.mean(d["seconds"]) / 2)
        row[name] = dict(sd_rep=float(E[:, k].std(ddof=1)), reps_in_B=float(n_rep),
                         rmse_at_B=float(E[:, k].std(ddof=1) / np.sqrt(n_rep)))
    row["lam_diff"] = dict(hybrid=hl["diff"][e], hybrid_se=hl["diff_se_crn"][e],
                           cond_rqmc_ref=ref["diff_lam05_minus_0"][e], cond_rqmc_ref_se=ref["diff_se_crn"][e],
                           cond_rqmc_small=cq["diff_lam05_minus_0"][e], cond_rqmc_small_se=cq["diff_se_crn"][e])
    out["rows"].append(row)
json.dump(out, open(f"{_ROOT}/out/r1_equal_time_v02.json", "w"), indent=1)
print(f"budget {B:.1f}s; cond RQMC per-rep per-lambda {out['cond_rqmc_sec_per_rep_per_lambda']:.2f}s")
for r in out["rows"]:
    h = r["hyb(3,3,2,2,2,2)"]; h3 = r["hyb(3,3,3,3,3,3)"]; h5 = r["hyb(3,3,3,5,5,5)"]
    print(f"{r['event']:28s} ref {r['ref']:.5f}±{r['ref_se']:.1e} | hyb332222 bias {h['bias']:+.1e} (z {h['bias_z']:+.1f}) rmse {h['rmse_at_B']:.1e}"
          f" | hyb333333 bias {h3['bias']:+.1e} | hyb333555 bias {h5['bias']:+.1e} | condRQMC rmse {r['cond_rqmc']['rmse_at_B']:.1e} | iid {r['cond_iid']['rmse_at_B']:.1e}"
          f" | dlam hyb {r['lam_diff']['hybrid']:+.4f}±{r['lam_diff']['hybrid_se']:.0e} ref {r['lam_diff']['cond_rqmc_ref']:+.4f}±{r['lam_diff']['cond_rqmc_ref_se']:.0e}")
