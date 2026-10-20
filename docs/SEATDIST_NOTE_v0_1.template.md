# Seat distributions for multi-party first-past-the-post elections by conditional exact recursion and factor integration: method and a numerical audit (Great Britain 2024 case)

Nari Han (Risk & Return Co., Ltd.; ORCID 0009-0008-7228-8740)

Draft v0.1, 2026-10-20. Not peer reviewed. Comments welcome.

## Abstract

Probabilistic seat forecasts for first-past-the-post (FPTP) elections are often computed by simulating the whole
election many times. If the model makes seat outcomes independent once a set of common factors (national polling
error, regional swings, seat-type factors) is fixed, the distribution of each party's seat total can instead be
computed exactly for every value of the factors, with the Poisson-binomial recursion familiar from credit portfolio
loss models, and only the factors need to be integrated numerically. This note describes that computation for a
multi-party model of Great Britain (631 seats, eight party categories, 66 factor dimensions) and, mainly, reports an
audit of its numerical error. In the scheme we used most, an outer randomised quasi-Monte Carlo (RQMC) integral over
nine factors wraps an inner Gauss–Hermite rule over regional factors. Its usual error bar, the spread over independent
scramblings, missed two errors of different kinds: a deterministic bias of the regional rule (P(Green ≥ 4 seats) off
by 3.2e-3 while its scramble standard error was 3.5e-5) and an under-resolved rare outer contribution (P(Other ≥ 5)).
We describe an audit that separates the error components (paired regional-rule comparisons at fixed outer points,
a conditional RQMC reference without regional quadrature, a regional-only RQMC bias estimate, and an independent
generative Monte Carlo). After the audit, a strengthened regional rule agrees with pooled references to within 1e-3 for
all 18 reported event probabilities (16 of 18 in a secondary, bias-corrected comparison), in an exploratory audit
with per-event nominal intervals. We do not claim a general
speed advantage over Monte Carlo, and the model is used here as a fixed computational case, not as a forecasting
method: it misjudged Labour's 2024 seat efficiency badly, which is studied separately.

## 1. Scope and claims

This is a computation note. Its subject is how to compute, and how to check, the seat distribution implied by a given
factor model. Three limits apply throughout.

- **Fixed computational case.** The Great Britain model (v0.2.1, final 2024 polls, no poll-bias shift, seat-type
  factors drawn from a prior built on 2010–2019 descriptive statistics) is held fixed. It is not presented as a good
  forecaster: with the actual 2024 vote shares, and seat-type factors set anywhere from zero to values derived from
  2024 itself, it gives Labour about 345 to 368 seats against 411 won. Those prior means
  are descriptive statistics, not estimates of generative coefficients.
- **No novelty claim for the recursion.** Exact seat or loss distributions under conditional independence are standard
  (Poisson-binomial computations, credit portfolio recursions, electoral-vote totals built from state-level probabilities). The contribution is
  the error audit for a multi-party, high-dimensional case and the documentation of where a common error bar fails.
- **No general efficiency claim.** Timings come from a shared two-core machine and are indicative only. We did not
  benchmark importance sampling, and we do not claim that deterministic integration beats Monte Carlo in general or
  for tail events.

The claims the note makes, and their support:

| # | claim | support |
|---|---|---|
| C1 | Under factor-conditional independence, each party's marginal seat distribution under the model is exact given the factors (recursion plus convolution across regions); only the factor integral is approximate. | §3; enumeration tests |
| C2 | In a hybrid scheme (outer RQMC, inner fixed quadrature), the outer scramble standard error does not contain the inner rule's bias. | Tables N2, N3, N5 (Green) |
| C3 | A few independent scramblings can understate a rare outer contribution; the estimate and its standard error were both unstable. | Table N4 (Other ≥ 5) |
| C4 | The error components can be separated by paired rule changes at fixed outer points, a conditional RQMC reference, a regional-only RQMC bias estimate and an independent generative Monte Carlo. In this case the procedure located both failures and led to a regional rule that agrees with pooled references within 1e-3 for 18/18 events (16/18 in the secondary, bias-corrected comparison; exploratory audit). | Tables N3–N7 |
| C5 | A target-event dimension-adaptive tensor search found a grid that was cheaper per run than our hand design and had a smaller maximum difference from Monte Carlo, but it was worse on 5 of 18 events and the search cost about 39 final runs; any gain from reusing the grid across scenarios is a hypothesis. | Table N9 |

## 2. Model class

Seats i = 1, …, n (here n = 631, Great Britain without the Speaker's seat); party categories k (Con, Lab, LD, Ref,
Grn, SNP in Scotland, PC in Wales, and Oth, the largest other candidate in the seat; the remaining minor vote is a
non-competing residual). The log vote index is

```
eta_ik = log b_ik + dbar_{B(i),k} + Delta_k + [A z]_k + [C_r(i) w_r(i)]_k + sum_j L_ikj mu_j + u_ik ,
```

where `b` is the baseline share (2019 notional result on 2024 boundaries, with the baseline adjustments of the model),
`dbar` the block (England, Scotland, Wales) swing that matches the polls, `Delta` an optional poll-bias shift (zero
here), `z ~ N(0, I_5)` national polling error mapped to vote-share covariance through `A`, `w_r` regional factors
(five per English region, six in Scotland and in Wales; 57 in all), `mu` four seat-type factors with seat-specific
loadings `L` (anti-Conservative consolidation, beyond-proportional Conservative change, Leave share on Reform, Muslim
share on Labour and independents), and `u` seat-level noise.

With Gumbel seat noise on a two-scale mixture (probability 1 − π with scale τ1, π with τ2; π = 0.176, τ1 = 0.114,
τ2 = 0.237, fitted to 2010–2019 residuals), the probability that party k wins seat i given all common factors
F = (z, w, mu) has a closed form:

```
p_ik(F) = (1 − π) softmax_k(eta_i / τ1) + π softmax_k(eta_i / τ2),     eta without u .
```

The one structural assumption the computation uses is **conditional independence**: given F, seat outcomes are
independent. All correlation between seats enters through the 66 common factor dimensions (5 national, 57 regional,
4 seat-type). Given z and mu, the regional factors of different regions are independent.

## 3. Conditional exact computation

**Recursion.** Given F, the seat total of party k is a sum of independent Bernoulli(p_ik(F)) variables, so its
distribution is Poisson-binomial. We build it seat by seat,

```
P_{j}(s) = P_{j−1}(s) (1 − p_jk) + P_{j−1}(s − 1) p_jk ,
```

which combines non-negative terms only (apart from forming 1 − p), so there is no cancellation and small tail
probabilities keep their relative precision unless they underflow (Chen & Liu 1997 survey the Poisson-binomial distribution; the same recursion is used for credit portfolio losses,
Andersen, Sidenius & Basu 2003; electoral-vote totals with margins of error have been estimated from state-level polls by Kaplan & Barnett 2003).

**Regions and convolution.** Given (z, mu), the seat distribution of each region is the integral over that region's
factors w_r of its conditional Poisson-binomial distribution; region distributions are then convolved. The national
distribution is the integral of the convolution over (z, mu).

**Events.** A majority (S_k ≥ 326) and similar events follow directly from party marginals. Majority events of
different parties are mutually exclusive, so P(hung parliament) = 1 − Σ_k P(S_k ≥ 326) exactly. The total of a group
of parties (for example Lab + LD) is computed by adding the seat-level win probabilities of the group before the
recursion; convolving party marginals would be wrong because the parties compete for the same seats. Events that
constrain two parties separately (for example Lab ≥ 300 and LD ≥ 60) need a two-dimensional recursion and are not
supported in this version.

**Truncation.** Seats whose conditional win probability is within ε = 1e-12 of 0 or 1 are treated as certain. We
report the resulting total-variation bound per factor node (here 3.1e-12 to 3.3e-12 for the hybrid runs of Table N3,
maximum over the evaluated nodes; the a-priori bound for the active-support kernel is n·ε = 6.3e-10 per party). It is
not a bound over the continuous factor space and not a bound on the total error.

**Sensitivities.** For a fixed F, ∂P(S_k ≥ m | F)/∂p_ik(F) = P(S_{k,−i} = m − 1 | F), the probability that seat i is
pivotal. Parameter sensitivities follow by integrating P(S_{k,−i} = m − 1 | F) ∂p_ik(F)/∂θ over F (for the mixture,
a weighted sum of the two softmax Jacobians). This identity is stated for completeness; it is not implemented in the
code described here.

**Implementation.** Python with Numba, parallel over factor nodes; exp(eta/τ) is factorised so that the inner loops
only multiply. An active-support variant of the recursion (skipping seats within ε of 0 or 1) gave the same party
PMFs as the dense recursion to 1.4e-14 and was about 2.4 times faster at the same factor points.

## 4. Integration schemes and error components

Four schemes compute the same quantity:

| scheme | national and seat-type factors (9 dims) | regional factors (57 dims) | seat outcomes | used as |
|---|---|---|---|---|
| S1 deterministic tensor (v0.1 engine) | rotated national axes: trapezoid on the leading axis, Gauss–Hermite on the rest (no seat-type factors in v0.1) | Gauss–Hermite tensor per region | exact recursion | main engine of v0.1 |
| S2 hybrid (v0.2.1 engine) | scrambled Sobol, 2^m points × R independent scramblings | Gauss–Hermite tensor per region, rule `reg_ns` (points on the regional axes ranked 1–6) | exact recursion | main engine of v0.2.1 |
| S3 conditional RQMC | scrambled Sobol over all 66 dimensions | (sampled with the rest) | exact recursion | reference |
| S4 generative Monte Carlo | sampled | sampled | sampled (Gumbel-max) | independent check of the implementation |

The error of S2 has three components with different behaviour:

1. **Outer randomisation error.** Estimated by the spread over independent scramblings (estimated SE). With few
   scramblings this estimate can be unstable, in particular when rare outer regions contribute (Table N4).
2. **Regional-rule bias.** Deterministic given the rule; it is common to all scramblings and therefore invisible in
   their spread (Tables N2, N5). Changing the rule at the same outer points gives paired differences that are a
   sensitivity measure, not an error bound: Gauss–Hermite rules of different orders are not nested, and the
   differences need not shrink monotonically (Table N3).
3. **Truncation.** Node-wise total-variation bound (§3), negligible here.

A reference from S3 is an estimate with its own sampling error, not the true value. We therefore audit with four
pieces: (a) paired regional-rule ladders at fixed outer points; (b) a conditional RQMC reference; (c) a direct
estimate of the regional-rule bias, comparing the hybrid core with regional-only conditional RQMC at the same outer
points, using independent outer scrambled nets so that the replicate spread is a valid standard error; and (d) an
independent generative Monte Carlo that shares no integration code with the engine.

## 5. Results for the Great Britain 2024 case

Eighteen event probabilities are reported throughout (a Labour majority, seat thresholds for each party, the actual
outcome or more/fewer, a hung parliament). They share the same distributions and are not independent: for example a Labour
majority and a hung parliament are close to complements (correlation −0.9997 across reference replicates), so the
number of events does not add up evidence.

### 5.1 Deterministic tensor quadrature (model v0.1) and its check by generative Monte Carlo

{{N1}}

The reference grid agrees with 3 × 10⁶ generative draws for 17 of 18 events (|z| ≤ 2.5). The exception, Green ≥ 4
(difference 2.1e-3), came from too few points on the regional axes that carry Green's convex response. A larger earlier
grid (24,705 national nodes, 769 s) had ordered the national axes by raw seat influence; the axis that moves the SNP ranked fourth
and received three points, giving errors of 1.4e-2 (SNP ≤ 9) and 1.6e-2 (Other ≥ 5). A grid ladder that refined only
some axes estimated its own error at 6e-4 and missed these; the independent Monte Carlo found them. Ordering the axes by
influence normalised per party fixed the allocation.

### 5.2 Hybrid scheme with the default regional rule

{{N2}}

Two events differ from the generative Monte Carlo by about three standard errors, while the hybrid's own scramble
standard error for Green is 3.5e-5. The next tables locate the two causes.

### 5.3 Regional rule at fixed outer points

{{N3}}

Raising the first two regional axes from 3 to 5 points moves every event by at most 1.1e-4. Raising axes 3–6 from 2 to
3 points moves Green by −3.4e-3 and Plaid Cymru by −5.3e-4; raising axes 4–6 further from 3 to 5 points moves Green
again by +6.9e-4 while the other 17 events move by at most 4.4e-5. The regional rule is therefore the source of the
Green discrepancy, the changes are not monotone, and the rule (3,3,3,3,3,3) is not yet converged for Green.

### 5.4 Outer resolution of a rare event

{{N4}}

Other ≥ 5 hardly moves under the regional rules tried (0.000339 in all of them). With four scramblings its estimate
was 0.000339 with a standard error of 1.1e-5; larger outer designs gave 0.00041 to 0.00050, consistent with each other
and with the reference 0.000453 ± 1.1e-5. Four scramblings had not revealed rare large outer contributions, so both the
estimate and its standard error were off. We do not attach normal-theory confidence statements to a handful of RQMC
replicates; Jain, Hickernell, Owen & Sorokin (2025) study intervals for RQMC that do not rely on them.

### 5.5 Comparison with the conditional RQMC reference

{{N5}}

Against the reference, the default hybrid is off by +3.2e-3 for Green (scramble SE 3.5e-5) and by +8.1e-4 for Plaid
Cymru (SE 1.3e-4). The regional rule (3,3,3,3,3,3) reduces the Green difference to −2.5e-4. Differences between two
poll-bias scenarios (λ = 0.5 versus 0) are smooth and precise in both schemes when common random numbers are used, so
we no longer claim that smooth scenario differences are specific to deterministic integration; the Green difference
still inherits the regional-rule bias (−0.1675 against −0.1645).

### 5.6 Green: which regional axes matter

{{N6}}

Axis 4 accounts for most of the Green change (+5.46e-4 from 3 to 9 points, converged at 9 points: 7 → 9 changes it by
1e-6); axis 5 adds +1.6e-4 to +1.7e-4; axis 6 exists only in Scotland and Wales (English regions have five regional factors) and does not affect Green or
Plaid Cymru (changes of order 1e-17 and 1e-14). The effects of axes 4 and 5 add up, so mixed refinements are not
needed here. Plaid Cymru is insensitive to these axes. Its difference in Table N5 came from the single outer run used
there (scramble SE 1.3e-4); it is absent from the audited comparison (Table N7a: −8.3e-5, interval containing zero;
Table N7b: regional-rule bias −8.5e-7 ± 2.5e-5).

### 5.7 The audited rule

We adopted the rule (3,3,3,9,5,3) and compared it with the pooled references under rules fixed before the additional
runs: three independent conditional RQMC references pooled with equal weight; the hybrid level taken from two
independent outer runs of 2¹² × 8 points; the paired increment from (3,3,3,3,3,3) to the audited rule added; per-event
nominal 95% t approximate intervals with Welch–Satterthwaite degrees of freedom; tolerance 1e-3 (the forecast reports event
probabilities to 0.1 percentage point), with 5e-4 shown for reference. The tolerance was chosen after point estimates had
been seen but before valid bias standard errors and the additional runs existed. This is an exploratory numerical
audit, not a pre-registered experiment.

{{N7a}}

The audited rule is within 1e-3 of the pooled references for all 18 events and within 5e-4 for 14. Three intervals
exclude zero (Con ≥ 121, Con ≤ 50, Grn ≥ 4), all below 3.5e-4, so small residual biases of that order are possible.
For Green, about +1.1e-4 of the +2.35e-4 difference is the estimated regional-rule bias:

{{N7b}}

Subtracting this estimate gives a secondary, bias-corrected table. Its intervals are wider because the bias estimate
adds variance, and two events (LD ≥ 50, Ref ≥ 50) can no longer be confirmed within 1e-3:

{{N7c}}

Two implementation details matter for reproducing this audit. First, the stored bias run seeded the outer net with
`default_rng([S, r])` and the inner scramblings with `default_rng([S, r, p, k])`; numpy's SeedSequence zero-pads short
entropy, so the first outer point's first inner stream coincided with the outer stream. The tables above drop the
first point of every net (seven points per net remain), which keeps the replicate means unbiased but is not a full
eight-point net; the code now separates the streams. Second, the code's default regional rule remains the fast
(3,3,2,2,2,2) so that stored runs reproduce; results in this note pass the audited rule explicitly.

### 5.8 Conditional RQMC versus conditional iid sampling

{{N8}}

At the same number of conditional evaluations, scrambled Sobol sampling of the 66 factors reduced the variance by
factors of 2.9 to 34 for 15 of the 18 events; the two smallest ratios (2.9 and 3.2) are close to the 97.5% point of
F(15, 15), about 2.86, so the gain is borderline for those two. For Lab ≤ 250, Ref ≥ 50 and Other ≥ 5 we observed no clear gain
(ratios 0.9–1.1). Ref ≥ 50 has probability about 0.33, so this is not only a rare-event effect; effective dimension and
the placement of factors in the Sobol sequence also matter. These are observed ratios from 16 replicates each.

## 6. Target-event dimension-adaptive tensor quadrature (exploratory)

For the v0.1 engine we tried to choose the grid automatically. The variables are the point counts on 11 allocation
axes (5 national axes and the regional axes ranked 1–6). From a coarse start, the search tries one more level on each
axis and adopts the axis whose refinement changes the 18 event probabilities most per second of run time; it stops
when the largest change falls below 4e-4 or a soft time budget is exhausted. The search uses no Monte Carlo; the Monte
Carlo column is added afterwards. This is a greedy variant of the dimension-adaptive idea of Gerstner & Griebel (2003)
on full tensor grids, with the target defined by the reported events.

{{N9}}

The final grid (19,250 nodes, 161 s per run, largest difference from Monte Carlo 1.2e-3) has a smaller maximum
difference than our hand design (22,575 nodes, 479 s, 2.1e-3), but it is not better event by event: it is worse on 5
of the 18 events, and three events differ from the 3 × 10⁶ Monte Carlo by more than 2.5 standard errors (Con ≥ 121
z = −4.1, Ref ≥ 50 z = +4.6, Grn ≥ 4 z = −3.4) against one for the hand design (Grn ≥ 4 z = +7.2). The trapezoid spacing on the leading national axis stayed at the coarsest
value (0.7) throughout, consistent in this case with the fast convergence of the trapezoid rule for smooth integrands
with Gaussian decay (Trefethen & Weideman 2014); we do not claim that result guarantees the accuracy here. Three
weaknesses remain: the search took 6,206 s, about 39 final runs; between 2,970 and 8,250 nodes the one-step criterion
fell while the error against Monte Carlo stayed near 3e-3, which suggests interactions between axes that a one-axis
look-ahead cannot see; and the stopping criterion was never reached. Whether reusing the grid across many scenarios
recovers the search cost is untested. A per-event stopping rule with absolute and relative tolerances, two-step
look-ahead or mixed refinements, and validation on scenarios not used in the search are the obvious next steps.

## 7. Practical recommendations

1. Report the integration scheme for each error component separately. For a hybrid scheme, the scramble standard
   error covers only the outer randomisation.
2. Audit the inner rule at fixed outer points (paired differences) and against a reference that does not use the inner
   rule. Treat rule differences as sensitivities, not bounds.
3. Use enough independent scramblings, and check rare events with larger outer designs; a small number of scramblings
   can understate both the estimate and its error.
4. When the bias of the inner rule is estimated by replicates, use independent outer nets and separate random streams
   for the outer and inner draws.
5. Keep an independent generative simulator that shares no integration code with the engine. In this project it
   flagged each of the allocation errors described above before a reference was available.
6. Fix pooling rules and tolerances before adding runs, and say when they were fixed.

## 8. Limitations

- One model, one election, one set of inputs. The 18 events are dependent.
- The references are themselves estimates. The final comparison is exploratory, with per-event nominal intervals and a
  tolerance chosen after seeing point estimates.
- Timings are from a shared two-core machine; no performance claims are made. Importance sampling (Glasserman & Li
  2005 for the credit analogue) was not benchmarked.
- Only party marginals and group totals are supported; joint constraints on two parties are not.
- The model's substantive shortcomings (seat efficiency in 2024, the meaning of the seat-type prior) are outside this
  note.

## 9. Reproducibility

Code, data and stored outputs are in the bundle that accompanies this note (folder `ukseat/`). Every table above is
generated from stored outputs by `scripts/seatdist_note_tables_v0_1.py`, and this document is assembled by
`scripts/build_seatdist_note_v0_1.py`.

| table | runs (scripts in `scripts/` → outputs in `out/`) |
|---|---|
| N1 | `ref_exact2_v0_1.py`, `ref_exact_v0_1.py`, `mc_compare_v0_1.py` → `ref_exact2_v0_1.json`, `ref_exact_v0_1.json`, `mc_compare_std_v0_1.json` |
| N2 | `v02_mc_check.py` → `v02_mc_check.json` |
| N3, N4 | `reg_convergence_v02.py`, `reg_convergence_v02b.py` → `reg_convergence_v02.json`, `reg_convergence_v02b.json` |
| N5 | `cond_mc_v02.py` (reference, seed 20261022), `hybrid_lam_crn_v02.py`; also reads `reg_convergence_v02.json` → `cond_rqmc_N8192_R32_s20261022.json`, `hybrid_lam_crn_v02.json` |
| N6 | `r1_followup_v02.py`, then `r1_green_axis4_v02.py`, then `r1_green_axis5_v02.py` (the last two add rows to the same file; running the first alone resets it) → `r1_followup_v02.json` |
| N7 | references: `cond_mc_v02.py` with seed 20261022 (run directly), 20261103 (called by `r1_green_resolve_v02.py`) and 20261104 (called by `r1_extra_runs_v02.py`); hybrid level: `r1_green_resolve_v02.py`, `r1_extra_runs_v02.py`; increment: `r1_followup_v02.json` (N6 chain); regional bias: `r1_fixed_outer_reps_v02.py`; tables: `r1_note_tables_v02.py` → `r1_note_tables_v02.json` |
| N8 | `cond_mc_v02.py` → `cond_rqmc_N2048_R16_s20261020.json`, `cond_iid_N2048_R16_s20261021.json` |
| N9 | `adaptive_quad_v0_1.py` → `adaptive_quad_v0_1.json` |

Two caveats for re-running. (i) The stored regional-bias run `out/r1_fixed_outer_reps_v02.json` was produced
before the random-stream fix (§5.7); the corrected script writes `out/r1_fixed_outer_reps_v02_seedfix.json` and does
not regenerate the stored file, whose re-aggregation (points 2–8 of each net) is what Tables N7 use. (ii) The stored
conditional runs used the dense recursion kernel (`UKSEAT_PB=dense`); the current default is the active-support kernel,
which agrees to 1.4e-14 in the party PMFs.

Requirements: Python 3.11 or later, numpy, scipy, pandas, numba. Unit tests: `python tests/test_engine_v0_1.py`
(enumeration of small elections against the recursion, truncation bound, quadrature moments, closed-form mixture
against Gumbel-max sampling).

Data: UK Parliament / House of Commons Library election results (psephology datasette, notional 2019 results on 2024
boundaries and 2024 results); British Election Study 2024 General Election results file (CC BY 4.0,
doi:10.48420/26340568); British Polling Council summary of final 2024 polls.

## 10. AI usage disclosure

Generative AI tools were used extensively. Claude (Anthropic; Claude Opus 5.5, October 2026) wrote most of the code,
ran the numerical experiments and drafted this note, working from the author's instructions. GPT 6.1 Sol (OpenAI)
reviewed the computations and drafts in several rounds and independently re-ran parts of the audit; changes made in
response are part of the audit record. The author chose the questions, directed and reviewed the work, decided which
results to report and how to word them, and takes full responsibility for the content.

## References

- Andersen, L., Sidenius, J. and Basu, S. (2003). All your hedges in one basket. *Risk*, November 2003, 67–72.
- Chen, S. X. and Liu, J. S. (1997). Statistical applications of the Poisson-binomial and conditional Bernoulli
  distributions. *Statistica Sinica* 7(4), 875–892.
- Gerstner, T. and Griebel, M. (2003). Dimension-adaptive tensor-product quadrature. *Computing* 71(1), 65–87.
  doi:10.1007/s00607-003-0015-5.
- Glasserman, P. and Li, J. (2005). Importance sampling for portfolio credit risk. *Management Science* 51(11),
  1643–1656. doi:10.1287/mnsc.1050.0415.
- Jain, A., Hickernell, F. J., Owen, A. B. and Sorokin, A. G. (2025). Empirical Bernstein and betting confidence
  intervals for randomized quasi-Monte Carlo. arXiv:2504.18677 (v2, 2026).
- Kaplan, E. H. and Barnett, A. (2003). A new approach to estimating the probability of winning the presidency.
  *Operations Research* 51(1), 32–40.
- Trefethen, L. N. and Weideman, J. A. C. (2014). The exponentially convergent trapezoidal rule. *SIAM Review* 56(3),
  385–458.
- British Election Study (2024). 2024 General Election results file, version 1.0. doi:10.48420/26340568.
- UK Parliament. psephology-datasette: election results data. https://github.com/ukparliament/psephology-datasette
- British Polling Council (2024). The performance of the polls in the 2024 general election.
  https://www.britishpollingcouncil.org/the-performance-of-the-polls-in-the-2024-general-election/
