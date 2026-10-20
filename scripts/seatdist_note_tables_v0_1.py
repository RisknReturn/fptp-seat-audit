"""Tables for the public computation note (SEATDIST_NOTE_v0_1.md), generated only from stored run outputs in out/.

Every number in the note's tables comes from this script. Sources:
  N1  v0.1 deterministic grid vs generative MC 3x10^6           ref_exact2_v0_1.json, ref_exact_v0_1.json, mc_compare_std_v0_1.json
  N2  v0.2.1 hybrid (default rule) vs generative MC              v02_mc_check.json
  N3  regional-rule ladder at fixed outer points (paired)        reg_convergence_v02.json, reg_convergence_v02b.json
  N4  outer resolution of a rare event                           reg_convergence_v02b.json
  N5  hybrid rules vs the conditional RQMC reference, lambda difference   reg_convergence_v02.json, cond_rqmc_N8192_R32_s20261022.json, hybrid_lam_crn_v02.json
  N6  Green: axis-by-axis regional audit                        r1_followup_v02.json
  N7  audited rule vs pooled references (final R1 tables)       r1_note_tables_v02.json
  N8  conditional RQMC vs conditional iid, same N               cond_rqmc_N2048_R16_s20261020.json, cond_iid_N2048_R16_s20261021.json
  N9  target-event dimension-adaptive search trace              adaptive_quad_v0_1.json
Writes out/seatdist_note_tables_v0_1.md.
"""
import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import json, numpy as np
OUT = f"{_ROOT}/out"
J = lambda f: json.load(open(f"{OUT}/{f}"))
SHORT = {"Lab>=326 (majority)": "Lab ≥ 326 (majority)", "Lab>=400": "Lab ≥ 400", "Lab>=450": "Lab ≥ 450",
         "Lab<=411 (actual or fewer)": "Lab ≤ 411", "Lab<=250": "Lab ≤ 250", "Con>=121 (actual or more)": "Con ≥ 121",
         "Con<100": "Con < 100", "Con<=50": "Con ≤ 50", "Con>=250": "Con ≥ 250", "LD>=72 (actual or more)": "LD ≥ 72",
         "LD>=50": "LD ≥ 50", "Ref>=5 (actual or more)": "Ref ≥ 5", "Ref>=50": "Ref ≥ 50", "Grn>=4 (actual or more)": "Grn ≥ 4",
         "SNP<=9 (actual or fewer)": "SNP ≤ 9", "PC>=4 (actual or more)": "PC ≥ 4", "Oth>=5 (actual or more)": "Oth ≥ 5",
         "Hung (no party>=326)": "Hung"}
EV = list(SHORT)
f4 = lambda x: f"{x:+.4f}".replace("-", "−")
L = []
p = lambda s="": L.append(s)


def e(x, sign=True, d=1):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "–"
    if x == 0:
        return "0"
    m, ex = f"{abs(x):.{d}e}".split("e")
    s = ("−" if x < 0 else "+") if sign else ("−" if x < 0 else "")
    return f"{s}{m}e{int(ex)}"


# ---------------- N1
ref, old, mc = J("ref_exact2_v0_1.json"), J("ref_exact_v0_1.json"), J("mc_compare_std_v0_1.json")
m6 = [r for r in mc["runs"] if r["N"] == 1_000_000]
p("### Table N1. Model v0.1: deterministic grid versus generative Monte Carlo (3 × 10⁶ draws)")
p()
p(f"Reference grid: national {ref['nodes']:,} nodes (trapezoid h = 0.35 on the leading rotated axis, Gauss–Hermite 7/5/5/3), regional GH (5,5,3,2,2,2); "
  f"{ref['sec']:.0f} s on 2 shared cores. Earlier reference grid: national axes ordered by raw seat influence ({old['nodes']:,} national nodes, {old['sec']:.0f} s).")
p()
p("| event | deterministic grid | MC 3 × 10⁶ | z | earlier grid | z |")
p("|---|---:|---:|---:|---:|---:|")
for k in EV:
    v = ref["ev"][k]; mv = float(np.mean([r["ev"][k] for r in m6])); se = np.sqrt(max(mv * (1 - mv), 1e-12) / 3e6); vo = old["ev"][k]
    p(f"| {SHORT[k]} | {v:.5f} | {mv:.5f} | {(v - mv) / se:+.1f} | {vo:.5f} | {(vo - mv) / se:+.1f} |".replace("| -", "| −"))
p()
# ---------------- N2
v2 = J("v02_mc_check.json")
p("### Table N2. Model v0.2.1, hybrid scheme with the default regional rule (3,3,2,2,2,2) versus generative MC (2 × 10⁵ draws)")
p()
p("| event | hybrid | scramble SE | MC | z |")
p("|---|---:|---:|---:|---:|")
for k in EV:
    x = v2[k]
    p(f"| {SHORT[k]} | {x['qmc']:.4f} | {e(x['qmc_se'], False)} | {x['mc']:.4f} | {x['z']:+.1f} |".replace("| -", "| −"))
p()
# ---------------- N3
rc, rb = J("reg_convergence_v02.json"), J("reg_convergence_v02b.json")
R = rc["rules"]; base = R["(3, 3, 2, 2, 2, 2)"]
cols = ["(3, 3, 3, 3, 3, 3)", "(5, 5, 2, 2, 2, 2)", "(5, 5, 3, 3, 3, 3)", "(5, 5, 5, 3, 3, 3)"]
pd_ = rb["paired_diffs"]["(3,3,3,5,5,5)-(3,3,3,3,3,3)"]
p("### Table N3. Regional rules at fixed outer points (2¹¹ × 4 scrambled Sobol, seed 0): paired differences")
p()
p("Rule = Gauss–Hermite points on the regional factor axes ranked 1–6. All runs share the outer points, so differences isolate the regional rule.")
p()
p("| event | default (3,3,2,2,2,2) | (3,3,3,3,3,3) − default | (5,5,2,2,2,2) − default | (5,5,3,3,3,3) − default | (5,5,5,3,3,3) − default | (3,3,3,5,5,5) − (3,3,3,3,3,3) [paired SE] | outer scramble SE (default) |")
p("|---|---:|---:|---:|---:|---:|---:|---:|")
for k in EV:
    row = [f"{base['mean'][k]:.6f}"] + [e(R[c]["diff_vs_default"][k]) for c in cols]
    row += [f"{e(pd_['mean'][k])} [{e(pd_['se'][k], False)}]", e(base["se"][k], False)]
    p(f"| {SHORT[k]} | " + " | ".join(row) + " |")
p()
p("Time (2 shared cores): " + ", ".join(f"{c.replace(' ', '')} {R[c]['sec']:.0f} s" for c in ["(3, 3, 2, 2, 2, 2)"] + cols) +
  f". Node-wise Poisson-binomial truncation bound (TV, max over nodes): {e(min(r['tv_max_node'] for r in R.values()), False)} to {e(max(r['tv_max_node'] for r in R.values()), False)}.")
p()
# ---------------- N4
p("### Table N4. Outer resolution of a rare event, P(Oth ≥ 5), default regional rule")
p()
p("| outer design | factor points | estimate | scramble SE | replicate range |")
p("|---|---:|---:|---:|---|")
p(f"| 2¹¹ × 4, seed 0 | 8,192 | {base['mean']['Oth>=5 (actual or more)']:.6f} | {e(base['se']['Oth>=5 (actual or more)'], False)} | – |")
for k, v in rb["outer"].items():
    c = v["cfg"]; n = (2 ** c["m_log2"]) * c["reps"]
    p(f"| 2^{c['m_log2']} × {c['reps']}, seed {c['seed']} | {n:,} | {v['mean']['Oth>=5 (actual or more)']:.6f} | {e(v['se']['Oth>=5 (actual or more)'], False)} | "
      f"{min(v['oth_reps']):.6f} – {max(v['oth_reps']):.6f} |")
cr = J("cond_rqmc_N8192_R32_s20261022.json")
p()
p(f"Conditional RQMC reference (no regional quadrature, {cr['N']:,} × {cr['reps']} scrambles): {cr['mean_0.0']['Oth>=5 (actual or more)']:.6f} ± {e(cr['se_0.0']['Oth>=5 (actual or more)'], False)} (estimated SE).")
p()
# ---------------- N5
hl = J("hybrid_lam_crn_v02.json")
p(f"### Table N5. Hybrid rules versus the conditional RQMC reference ({cr['N']:,} points × {cr['reps']} scrambles, {cr['dims']} factor dimensions, all sampled)")
p()
p("Single runs; a difference contains the sampling error of both the hybrid run and the reference. It is not an expected error or an RMSE.")
p()
p("| event | reference ± est. SE | default − ref | default scramble SE | (3,3,3,3,3,3) − ref | λ-difference (0.5 − 0): hybrid ± SE | λ-difference: reference ± SE |")
p("|---|---|---:|---:|---:|---:|---:|")
for k in EV:
    r0, s0 = cr["mean_0.0"][k], cr["se_0.0"][k]
    p(f"| {SHORT[k]} | {r0:.5f} ± {e(s0, False)} | {e(base['mean'][k] - r0)} | {e(base['se'][k], False)} | {e(R['(3, 3, 3, 3, 3, 3)']['mean'][k] - r0)} | "
      f"{f4(hl['diff'][k])} ± {e(hl['diff_se_crn'][k], False, 0)} | {f4(cr['diff_lam05_minus_0'][k])} ± {e(cr['diff_se_crn'][k], False, 0)} |")
p()
p(f"Measured times (2 shared cores, indicative only): default hybrid {base['sec']:.0f} s for one λ in the rule ladder, {hl['0.0']['sec']:.0f} s and {hl['0.5']['sec']:.0f} s in the two λ runs of the CRN comparison; (3,3,3,3,3,3) {R['(3, 3, 3, 3, 3, 3)']['sec']:.0f} s; reference {sum(cr['seconds']):.0f} s (dense recursion kernel).")
p()
# ---------------- N6
fu = J("r1_followup_v02.json")
p("### Table N6. Green: axis-by-axis regional audit at the same outer points (baseline (3,3,3,3,3,3))")
p()
p("| rule | Grn ≥ 4 − baseline [paired SE] | PC ≥ 4 − baseline | Grn ≥ 4 | time |")
p("|---|---:|---:|---:|---:|")
for rule, g in fu["green_axes"].items():
    k, kp = "Grn>=4 (actual or more)", "PC>=4 (actual or more)"
    p(f"| {rule.replace(' ', '')} | {e(g['diff_vs_333333'][k], d=3)} [{e(g['diff_se'][k], False)}] | {e(g['diff_vs_333333'][kp])} | {g['mean'][k]:.6f} | {g['sec']:.0f} s |")
p()
# ---------------- N7
nt = J("r1_note_tables_v02.json")
for key, title in (("primary", "Table N7a. Audited rule (3,3,3,9,5,3) minus the pooled references (three conditional RQMC references, 96 scrambles)"),
                   ("corrected", "Table N7c. As N7a after subtracting the estimated regional-rule bias (secondary)")):
    p(f"### {title}")
    p()
    p("Per-event nominal 95% t approximate intervals (Welch–Satterthwaite df); not simultaneous over the 18 events.")
    p()
    p("| event | pooled reference ± SE | difference | 95% interval (df) | within 1e-3 | within 5e-4 |")
    p("|---|---|---:|---:|:-:|:-:|")
    for k in EV:
        r = nt["events"][k]; x = r[key]
        p(f"| {SHORT[k]} | {r['ref_pooled']:.6f} ± {e(r['ref_se'], False)} | {e(x['diff_vs_ref'])} | [{e(x['ci95'][0])}, {e(x['ci95'][1])}] ({x['df']:.0f}) | "
          f"{'yes' if x['within']['0.001'] else '**no**'} | {'yes' if x['within']['0.0005'] else 'no'} |")
    n = nt["summary"][key + "_n_within"]
    p()
    p(f"Within 1e-3: {n['0.001']}/18; within 5e-4: {n['0.0005']}/18.")
    p()
p("### Table N7b. Regional-rule bias of the audited rule: hybrid minus regional-only conditional RQMC at the same outer points")
p()
p("16 independent outer scrambled Sobol nets; mean of points 2–8 of each net (the first point is dropped because its inner random stream coincided with the outer one in the stored run); SE from the 16 replicate means, t₁₅.")
p()
p("| event | bias ± SE | 95% interval |")
p("|---|---:|---:|")
from scipy.stats import t as _t
q = _t.ppf(0.975, 15)
for k in EV:
    r = nt["events"][k]; b, s = r["regional_bias_7pt"], r["regional_bias_7pt_se"]
    p(f"| {SHORT[k]} | {e(b)} ± {e(s, False)} | [{e(b - q * s)}, {e(b + q * s)}] |")
p()
# ---------------- N8
c1, c2 = J("cond_rqmc_N2048_R16_s20261020.json"), J("cond_iid_N2048_R16_s20261021.json")
Rr, Ri = np.array(c1["reps_events"]["0.0"]), np.array(c2["reps_events"]["0.0"])
ratio = Ri.var(0, ddof=1) / Rr.var(0, ddof=1)
p(f"### Table N8. Conditional sampling: observed variance ratio iid / RQMC at the same N ({c1['N']:,} points, {c1['reps']} replicates each)")
p()
p("| event | ratio | event | ratio | event | ratio |")
p("|---|---:|---|---:|---|---:|")
for i in range(0, 18, 3):
    p("| " + " | ".join(f"{SHORT[c1['events'][j]]} | {ratio[j]:.1f}" for j in range(i, i + 3)) + " |")
p()
# ---------------- N9
aq = J("adaptive_quad_v0_1.json")
p("### Table N9. Target-event dimension-adaptive tensor quadrature (model v0.1): search trace")
p()
p("Greedy rule: refine the axis with the largest change of the 18 event probabilities per second. MC is not used for selection or stopping; the MC column is shown afterwards.")
p()
p("| step | grid (national; regional) | nodes | time per run | max \\|grid − MC\\| | cumulative search time |")
p("|---:|---|---:|---:|---:|---:|")
import ast


def spec_txt(sp):
    nat, reg = ast.literal_eval(sp)
    tr = [x for x in nat if x[0] == "trap"][0]
    gh = "/".join(str(x[1]) for x in nat if x[0] == "gh")
    return f"trapezoid h = {tr[1]}, GH {gh}; regional ({','.join(map(str, reg))})"


for i, r in enumerate(aq["log"]):
    p(f"| {i} | {spec_txt(r['spec'])} | {r['nodes']:,} | {r['sec']:.1f} s | {r['err_vs_mc']:.4f} | {r['elapsed']:,.0f} s |")
f = aq["final"]
mcv = aq["mc"]
err_final = max(abs(f["ev"][k] - mcv[k]) for k in mcv)
p(f"| final | {spec_txt(f['spec'])} | {f['nodes']:,} | {f['sec']:.0f} s | {err_final:.4f} | – |")
err_hand = max(abs(ref["ev"][k] - float(np.mean([r["ev"][k] for r in m6]))) for k in ref["ev"])
p(f"| hand design (N1) | trapezoid h = 0.35, GH 7/5/5/3; regional (5,5,3,2,2,2) | {ref['nodes']:,} | {ref['sec']:.0f} s | {err_hand:.4f} | – |")
p()
open(f"{OUT}/seatdist_note_tables_v0_1.md", "w").write("\n".join(L) + "\n")
print("\n".join(L[:12])); print("...", len(L), "lines")
