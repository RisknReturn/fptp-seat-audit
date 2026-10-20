"""Plain Monte Carlo vs exact: accuracy on events (incl. deep tails) and wall-clock time,
plus noise in scenario differences (common random numbers)."""
import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import sys, json, time, numpy as np
sys.path.insert(0, _ROOT)
from ukseat.setup2024 import make_model
from ukseat.engine import monte_carlo
from ukseat.report import events_mc

m = make_model()
monte_carlo(m, 0.0, 200)                     # compile
plan = [(1_000, 30), (10_000, 30), (100_000, 10), (1_000_000, 3)]
if len(sys.argv) > 1 and sys.argv[1] == "big":
    plan = [(10_000_000, 1)]
out = []
OFFSET = {1_000: 10**9, 10_000: 2 * 10**9, 100_000: 3 * 10**9, 1_000_000: 4 * 10**9, 10_000_000: 5 * 10**8}
for N, reps in plan:
    for rep in range(reps):
        S, sec = monte_carlo(m, 0.0, N, seed=OFFSET[N] + rep * N)      # disjoint seed blocks
        ev = events_mc(S)
        out.append(dict(N=N, rep=rep, sec=sec, ev=ev, mean=S.mean(0).tolist(), sd=S.std(0).tolist()))
        print(f"N={N:>9d} rep={rep:2d} sec={sec:7.2f} Lab>=326={ev['Lab>=326 (majority)']:.5f} Lab<=250={ev['Lab<=250']:.2e}", flush=True)
# scenario differences with common random numbers: lambda 0.5 vs 0.25
diffs = []
for rep in range(10):
    seed = 6 * 10**8 + rep * 100_000                 # disjoint across replicates; same within a pair (CRN)
    a, _ = monte_carlo(m, 0.25, 100_000, seed=seed); b, _ = monte_carlo(m, 0.5, 100_000, seed=seed)
    ea, eb = events_mc(a), events_mc(b)
    diffs.append({k: eb[k] - ea[k] for k in ea})
tag = "big" if len(sys.argv) > 1 else "std"
json.dump(dict(runs=out, crn_diffs_N1e5=diffs), open(f"{_ROOT}/out/mc_compare_{tag}_v0_1.json", "w"))
print("done")
