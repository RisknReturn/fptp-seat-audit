"""Regenerate the R1 note's comparison table from raw outputs only (reply to the R1 follow-up review, item 4).

Inputs (raw):
  out/cond_rqmc_N8192_R32_s20261022.json, out/cond_rqmc_N8192_R32_s20261103.json   two conditional RQMC references
  out/r1_green_resolve_v02.json  ["outer_seed1_333333"]   hybrid rule (3,3,3,3,3,3), outer 2^12 x 8 scrambles, seed 1
  out/reg_convergence_v02b.json  ["paired"]                hybrid (3,3,3,3,3,3), outer 2^11 x 4, seed 0 (secondary)
  out/r1_followup_v02.json       ["green_axes"]            paired increment (3,3,3,3,3,3) -> (3,3,3,9,5,3), seed 0
  out/r1_fixed_outer_reps_v02.json                         regional-rule bias of (3,3,3,9,5,3), 16 independent outer nets
      (re-aggregated from points[1:] of each replicate, i.e. 7-point means: the first outer point's inner RNG stream
       coincided with the outer stream in the stored run; this differs from a full 8-point Sobol net audit)
Rules fixed before the tables were produced (2026-10-20):
  * reference = all scramble replicates of all available references (two, or three once
    out/cond_rqmc_N8192_R32_s20261104.json exists), equal weight (SE from the replicate values);
  * hybrid level = the larger fixed-budget run(s) (2^12 x 8: seed 1, and seed 2 once out/r1_extra_runs_v02.json has it,
    equal weight) for EVERY event; the seed 0 run is shown only as a secondary check. No inverse-variance pooling of runs
    (its SEs were unstable for Oth and PC);
  * PRIMARY estimate for the audited rule = level(3,3,3,3,3,3) + paired increment, compared directly with the reference;
    SECONDARY (bias-corrected) estimate additionally subtracts the estimated regional bias (both tables are reported);
  * interval: per-event nominal 95% t approximate interval, Welch-Satterthwaite df (not simultaneous over events);
  * this is an exploratory numerical audit in which the pooling rules were fixed before the additional runs; it is not
    a pre-registered experiment;
  * tolerance: eps = 1e-3 (absolute probability; the note reports event probabilities to 0.1 percentage point),
    also shown against 5e-4. eps was set after the point estimates were seen but before the valid bias SEs existed.
  * z^2 means are descriptive only (events share the same PMFs; e.g. Lab majority vs Hung are almost complements).
Writes out/r1_note_tables_v02.json and out/r1_note_tables_v02.md."""
import os as _os, json, numpy as np
from scipy.stats import t as student_t
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
O = lambda f: json.load(open(f"{_ROOT}/out/{f}"))
r1, r2 = O("cond_rqmc_N8192_R32_s20261022.json"), O("cond_rqmc_N8192_R32_s20261103.json")
EV = r1["events"]; assert EV == r2["events"]
REFS = [r1, r2]
if _os.path.exists(f"{_ROOT}/out/cond_rqmc_N8192_R32_s20261104.json"):
    REFS.append(O("cond_rqmc_N8192_R32_s20261104.json")); assert REFS[-1]["events"] == EV
s1 = O("r1_green_resolve_v02.json")["outer_seed1_333333"]
s2 = O("r1_extra_runs_v02.json").get("outer_seed2_333333") if _os.path.exists(f"{_ROOT}/out/r1_extra_runs_v02.json") else None
rg = O("reg_convergence_v02b.json")["paired"]["(3, 3, 3, 3, 3, 3)"]
inc = O("r1_followup_v02.json")["green_axes"]["(3, 3, 3, 9, 5, 3)"]
fb = O("r1_fixed_outer_reps_v02.json"); assert fb["events"] == EV and fb["rule"] == [3, 3, 3, 9, 5, 3]
EPS = [1e-3, 5e-4]
A1 = np.array(r1["reps_events"]["0.0"]); A2 = np.array(r2["reps_events"]["0.0"])
AP = np.vstack([np.array(x["reps_events"]["0.0"]) for x in REFS])
S0 = np.array(rg["reps"])
# Regional bias from the stored replicate run, re-aggregated WITHOUT each replicate's first outer point (its inner RNG
# stream coincided with the outer stream; see the SEED FIX note in r1_fixed_outer_reps_v02.py): 7-point means per replicate.
B = np.array([np.mean([np.array(pt["hybrid"]) - np.array(pt["mc"]) for pt in x["points"][1:]], axis=0) for x in fb["replicates"]])


def interval(comps, d):
    v = sum(s ** 2 for s, _ in comps)
    df = v ** 2 / sum((s ** 4) / f for s, f in comps if s > 0) if v > 0 else np.inf
    half = student_t.ppf(0.975, df) * np.sqrt(v)
    return float(np.sqrt(v)), float(df), [d - half, d + half], {str(e): bool(abs(d) + half <= e) for e in EPS}


rows, out = [], {}
for i, k in enumerate(EV):
    ref, ref_se, ref_df = AP[:, i].mean(), AP[:, i].std(ddof=1) / np.sqrt(len(AP)), len(AP) - 1
    z12 = (A2[:, i].mean() - A1[:, i].mean()) / np.hypot(A1[:, i].std(ddof=1) / np.sqrt(len(A1)), A2[:, i].std(ddof=1) / np.sqrt(len(A2)))
    if s2 is None:
        h, h_se, h_df = s1["mean"][k], s1["se"][k], 8 - 1
    else:
        h = 0.5 * (s1["mean"][k] + s2["mean"][k]); h_se = 0.5 * np.hypot(s1["se"][k], s2["se"][k]); h_df = 14
    h0, h0_se = S0[:, i].mean(), S0[:, i].std(ddof=1) / np.sqrt(len(S0))
    dinc, dinc_se, dinc_df = inc["diff_vs_333333"][k], inc["diff_se"][k], 4 - 1
    b, b_se, b_df = B[:, i].mean(), B[:, i].std(ddof=1) / np.sqrt(len(B)), len(B) - 1
    # PRIMARY: the audited rule's own estimate against the reference (no regional-bias subtraction)
    d_u = h + dinc - ref
    se_u, df_u, ci_u, w_u = interval([(h_se, h_df), (dinc_se, dinc_df), (ref_se, ref_df)], d_u)
    # SECONDARY: estimate after subtracting the estimated regional bias (a separate, bias-corrected estimator)
    d_c = h + dinc - b - ref
    se_c, df_c, ci_c, w_c = interval([(h_se, h_df), (dinc_se, dinc_df), (b_se, b_df), (ref_se, ref_df)], d_c)
    rec = dict(n_references=len(REFS), n_hybrid_runs=1 if s2 is None else 2,
               ref_pooled=ref, ref_se=ref_se, z_ref1_vs_ref2=z12, hybrid333_seed1=s1["mean"][k],
               hybrid333_seed2=None if s2 is None else s2["mean"][k], hybrid333_level=h, hybrid333_level_se=h_se,
               hybrid333_seed0=h0, hybrid333_seed0_se=h0_se, z_seed1_vs_seed0=(s1["mean"][k] - h0) / np.hypot(s1["se"][k], h0_se),
               increment_to_audited=dinc, increment_se=dinc_se, regional_bias_7pt=b, regional_bias_7pt_se=b_se,
               primary=dict(diff_vs_ref=d_u, se=se_u, df=df_u, ci95=ci_u, within=w_u),
               corrected=dict(diff_vs_ref=d_c, se=se_c, df=df_c, ci95=ci_c, within=w_c))
    out[k] = rec
    rows.append(rec)
summary = dict(n_events=len(EV), interval_note="per-event nominal 95% t approximate intervals (Welch-Satterthwaite df); not simultaneous",
               primary_n_within={str(e): int(sum(r["primary"]["within"][str(e)] for r in rows)) for e in EPS},
               corrected_n_within={str(e): int(sum(r["corrected"]["within"][str(e)] for r in rows)) for e in EPS},
               corr_LabMaj_Hung_in_reference_reps=float(np.corrcoef(AP[:, EV.index("Lab>=326 (majority)")],
                                                                    AP[:, EV.index("Hung (no party>=326)")])[0, 1]))
json.dump(dict(rules=__doc__, events=out, summary=summary), open(f"{_ROOT}/out/r1_note_tables_v02.json", "w"), indent=1)
L = []
for key, label in (("primary", "PRIMARY: audited rule (3,3,3,9,5,3) estimate - pooled reference (no regional-bias subtraction)"),
                   ("corrected", "SECONDARY: after subtracting the estimated regional bias (7-point replicate means, first outer point dropped)")):
    L += [f"### {label}", "", "per-event nominal 95% t approximate intervals; not simultaneous intervals over the 18 events", "",
          "| event | pooled reference ± SE | difference | 95% interval (df) | within 1e-3 | within 5e-4 |",
          "|---|---:|---:|---:|:-:|:-:|"]
    for k in EV:
        r = out[k]; x = r[key]
        L.append(f"| {k} | {r['ref_pooled']:.6f} ± {r['ref_se']:.1e} | {x['diff_vs_ref']:+.1e} | "
                 f"[{x['ci95'][0]:+.1e}, {x['ci95'][1]:+.1e}] ({x['df']:.0f}) | "
                 f"{'yes' if x['within']['0.001'] else 'no'} | {'yes' if x['within']['0.0005'] else 'no'} |")
    n = summary[key + "_n_within"]
    L += ["", f"within 1e-3: {n['0.001']} / {len(EV)}; within 5e-4: {n['0.0005']} / {len(EV)}", ""]
L += ["regional bias (7-point replicate means, 16 replicates): " + "; ".join(f"{k.split(' ')[0]} {out[k]['regional_bias_7pt']:+.1e} ± {out[k]['regional_bias_7pt_se']:.1e}" for k in EV),
      f"corr(Lab majority, Hung) over reference replicates {summary['corr_LabMaj_Hung_in_reference_reps']:.4f}"]
open(f"{_ROOT}/out/r1_note_tables_v02.md", "w").write("\n".join(L) + "\n")
print("\n".join(L))
