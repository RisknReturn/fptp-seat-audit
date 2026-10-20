"""Additional table for the computation note v0.2.2 (review of 2026-10-21), generated only from stored outputs in out/.
  N7d  audited rule run directly at fresh outer points (seed 3000), vs the pooled reference of N7a
       sources: audited_fresh_seed_v0_1.json, r1_note_tables_v02.json
Tables M1 and N1-N9 are unchanged and still come from scripts/seatdist_note_tables_v0_2.py.
Writes out/seatdist_note_tables_v0_2_2.md."""
import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import sys, json
sys.path.insert(0, f"{_ROOT}/scripts")
OUT = f"{_ROOT}/out"
J = lambda f: json.load(open(f"{OUT}/{f}"))
# shared labels and number format of the v0.2 tables (import without running that script's table code)
src = open(f"{_ROOT}/scripts/seatdist_note_tables_v0_2.py").read()
ns = {"__file__": f"{_ROOT}/scripts/seatdist_note_tables_v0_2.py"}
exec(src[: src.index("# ---------------- M1")], ns)
SHORT, EV, e = ns["SHORT"], ns["EV"], ns["e"]

fr, nt = J("audited_fresh_seed_v0_1.json"), J("r1_note_tables_v02.json")
d = fr["design"]
L = []
p = lambda s="": L.append(s)
p("### Table N7d. Audited rule (3,3,3,9,5,3) run directly at fresh outer points, minus the pooled reference of Table N7a")
p()
p(f"Hybrid: rule fixed in advance, 2¹² × {d['reps']} scrambled Sobol, seed {d['seed']} (not used before), λ = 0, no paired increment; "
  "reference and interval method as in Table N7a (per-event nominal 95% t approximate intervals, Welch–Satterthwaite df; "
  "not simultaneous). No pass/fail rule was attached to this run.")
p()
p("| event | pooled reference ± SE | direct run ± scramble SE | difference | 95% interval (df) | N7a difference | within 1e-3 | within 5e-4 |")
p("|---|---|---|---:|---:|---:|:-:|:-:|")
for k in EV:
    r, a = fr["events"][k], nt["events"][k]["primary"]
    p(f"| {SHORT[k]} | {r['ref']:.6f} ± {e(r['ref_se'], False)} | {r['hybrid']:.6f} ± {e(r['hybrid_se'], False)} | {e(r['diff'])} | "
      f"[{e(r['ci95'][0])}, {e(r['ci95'][1])}] ({r['df']:.0f}) | {e(a['diff_vs_ref'])} | "
      f"{'yes' if r['within_1e3'] else '**no**'} | {'yes' if r['within_5e4'] else 'no'} |")
p()
p(f"Within 1e-3: {fr['n_within_1e3']}/18; within 5e-4: {fr['n_within_5e4']}/18. Run time {fr['sec'] / 60:.0f} min on two cores; "
  f"largest per-node truncation bound {e(fr['tv_max_node'], False)}.")
p()
open(f"{OUT}/seatdist_note_tables_v0_2_2.md", "w").write("\n".join(L))
print("\n".join(L))
