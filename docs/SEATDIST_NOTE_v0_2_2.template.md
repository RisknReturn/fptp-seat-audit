# A numerical audit of factor-based FPTP seat distributions: a Great Britain 2024 case study

Nari Han (Risk & Return Co., Ltd.; ORCID 0009-0008-7228-8740)

Draft v0.2.2, 2026-10-21. Not peer reviewed. Comments welcome. Table numbers N1–N9 are kept from v0.1 for cross-reference.
v0.2.2 adds a direct run of the audited rule at fresh outer points (Table N7d), a note on relative accuracy for the two
rarest events (Section 5.4) and two code fixes that change no reported number (Appendix F).

## Abstract

Probabilistic seat forecasts for first-past-the-post (FPTP) elections are often computed by simulating the whole
election many times. If the model makes seat outcomes independent once a set of common factors (national polling
error, regional swings, seat-type factors) is fixed, the distribution of each party's seat total can instead be
computed for every value of the factors with the Poisson-binomial recursion familiar from credit portfolio loss
models. Conditional on the common factors, the recursion and the convolution across regions are algebraically exact;
the implementation additionally incurs floating-point error and a bounded support-truncation error, and the integral
over the factors is the dominant approximation audited here. This note reports such an audit for a fixed multi-party
model of Great Britain (631 seats, eight party categories, 66 factor dimensions). In the scheme we used most, an outer
randomised quasi-Monte Carlo (RQMC) integral over nine factors wraps an inner Gauss–Hermite rule over regional
factors. Its usual error bar, the spread over independent scramblings, missed two errors of different kinds: a
deterministic bias of the regional rule (P(Green ≥ 4 seats) off by about 3e-3 while its scramble standard error was
3.5e-5) and an under-resolved rare outer contribution (P(Other ≥ 5)). We separate the error components with paired
regional-rule comparisons at fixed outer points, a conditional RQMC reference without regional quadrature, a
regional-only RQMC bias estimate and an independent generative Monte Carlo. With a strengthened regional rule, for all
18 reported events the per-event nominal 95% intervals for the estimated hybrid–reference differences lie within
[−1e-3, 1e-3] (16 of 18 in a secondary, bias-corrected comparison; 18 of 18 again in a direct run of the fixed rule at
new outer points). This is an exploratory comparison, not a simultaneous or deterministic error guarantee. We do not
claim a general speed advantage over Monte Carlo, and the model is used as a fixed computational case, not as a
forecasting method: it misjudged Labour's 2024 seat efficiency badly, which is studied separately.

## 1. Scope and claims

This is a computation note. Its subject is how to compute, and how to check, the seat distribution implied by a given
factor model. Three limits apply throughout.

- **Fixed computational case.** The main audit (Sections 4–6) concerns one fixed model: v0.2.1, final 2024 polls, no
  poll-bias shift, seat-type factors drawn from the prior of Section 2. Results for the earlier model v0.1 (no
  seat-type factors) are a separate, earlier computation, reported in Appendix A; probability differences between the
  two models are not differences of integration accuracy. The model is not presented as a good forecaster: with the
  actual 2024 vote shares, and seat-type factors set anywhere from zero to values derived from 2024 itself, it gives
  Labour about 345 to 368 seats against 411 won (the upper end uses factors derived from 2024 results). The prior
  means of the seat-type factors are descriptive statistics, not estimates of generative coefficients.
- **No novelty claim for the recursion or for the general principles.** Exact seat or loss distributions under
  conditional independence are standard (Poisson-binomial computations, credit portfolio recursions, electoral-vote
  totals built from state-level probabilities), and so are the facts that an outer standard error does not see an
  inner rule's bias and that a few randomisations can mislead. The contribution is a reproducible case: two failures
  of a common error bar in a high-dimensional multi-party model, each located and checked by a different computation
  path.
- **No general efficiency claim.** Timings come from a shared two-core machine and are indicative only. We did not
  benchmark importance sampling, and we do not claim that deterministic integration beats Monte Carlo in general or
  for tail events.

The claims the note makes, and their support:

| # | claim | support |
|---|---|---|
| C1 | Under factor-conditional independence, the Poisson-binomial recursion and the convolution across regions give each party's seat distribution given the factors algebraically exactly. The implementation adds floating-point error and a bounded support-truncation error; the factor integral is the dominant approximation. | §3; enumeration tests |
| C2 | In a hybrid scheme (outer RQMC, inner fixed quadrature), the outer scramble standard error does not estimate the inner rule's mean bias. | Table M1; Tables N2, N3, N5 (Green) |
| C3 | With a few independent scramblings, the estimate of a rare outer contribution and its standard error were both unstable in this case. | Table M1; Table N4 (Other ≥ 5) |
| C4 | The error components can be separated by paired rule changes at fixed outer points, a conditional RQMC reference, a regional-only RQMC bias estimate and an independent generative Monte Carlo. In this case the procedure located both failures. For the strengthened rule, the per-event nominal 95% intervals for the hybrid–reference differences lie within [−1e-3, 1e-3] for 18/18 events (16/18 in the secondary, bias-corrected comparison; 18/18 in a direct run of the fixed rule at fresh outer points). This is an exploratory per-event comparison, not a simultaneous error guarantee. | Tables N7a, N7d; Appendix C |
| C5 | For the earlier model v0.1, a target-event dimension-adaptive tensor search found a grid that was cheaper per run than our hand design and had a smaller maximum difference from Monte Carlo, but it was worse on 5 of 18 events and the search cost about 39 final runs; any gain from reusing the grid across scenarios is a hypothesis. | Appendix A (Table N9) |

## 2. Model class

Seats i = 1, …, n (here n = 631: the 632 Great Britain seats without the Speaker's seat; the 18 Northern Ireland
seats are outside the model). Party categories k: Con, Lab, LD, Ref, Grn, SNP (Scotland), PC (Wales), and Oth, the
largest other candidate in the seat (independents and minor parties together form this one category; the remaining
minor vote is a non-competing residual). The log vote index is

```
eta_ik = log b_ik + dbar_{B(i),k} + Delta_k + [A z]_k + [C_r(i) w_r(i)]_k + sum_j L_ikj mu_j + u_ik ,
```

where `b` is the baseline share (2019 notional result on 2024 boundaries, with the baseline adjustments of the model),
`dbar` the block (England, Scotland, Wales) swing that matches the polls, `Delta` an optional poll-bias shift (zero
here), `z ~ N(0, I_5)` national polling error mapped to vote-share covariance through `A`, `w_r ~ N(0, I)` regional
factors (five per English region, six in Scotland and in Wales; 57 in all), `mu` four seat-type factors with
seat-specific loadings `L`, and `u` seat-level noise.

**Seat-type factors in the audited model (v0.2.1).** The four factors are independent normal variables, independent of
z and w:

| factor | loading (before centring) | prior mean m_j | prior sd s_j |
|---|---|---:|---:|
| mu_T anti-Conservative consolidation | +1 on the stronger of Lab/LD and −1 on the other, in seats in England and Wales won by the Conservatives in the 2019 notional result | −0.0149 | 0.0604 |
| mu_C beyond-proportional Conservative change | Conservative baseline share minus its regional mean | −0.1200 | 0.3926 |
| mu_R Leave share on Reform | (Leave % − 52)/10 on Reform | 0.0480 | 0.0689 |
| mu_L Muslim share on Labour and independents | (Muslim % − mean)/10 on Labour; −κ·Muslim %/10 on Oth in seats with a Workers Party candidate or another independent or minor contender who is not a sitting MP standing as an independent, κ = 4 | 0.0174 | 0.0473 |

m_j and s_j are the mean and the sample standard deviation (divisor 2) of the factor values measured for the three
transitions 2010–15, 2015–17 and 2017–19 (`data/calib_mu_v02.json`). The loadings are centred per party with vote
weights and then adjusted so that the block vote shares do not move to first order at mu = 0 (`ukseat/model2.py`,
`neutral=True`).

**Seat noise.** For each seat one scale is drawn, τ2 with probability π and τ1 otherwise (π = 0.176, τ1 = 0.114,
τ2 = 0.237, fitted to 2010–2019 residuals), and is shared by all candidates in the seat; given the scale, u_ik = τ g_ik
with independent standard Gumbel g_ik. The probability that party k wins seat i given all common factors F = (z, w, mu)
then has a closed form:

```
p_ik(F) = (1 − π) softmax_k(eta_i / τ1) + π softmax_k(eta_i / τ2),     eta without u .
```

(A mixture with a separate scale per candidate would not give this formula. The generative Monte Carlo of Section 4
draws one scale per seat.)

The one structural assumption the computation uses is **conditional independence**: given F, seat outcomes are
independent. All correlation between seats enters through the 66 common factor dimensions (5 national, 57 regional,
4 seat-type). Given z and mu, the regional factors of different regions are independent.

## 3. Conditional computation

**Recursion.** Given F, the seat total of party k is a sum of independent Bernoulli(p_ik(F)) variables, so its
distribution is Poisson-binomial. We build it seat by seat,

```
P_{j}(s) = P_{j−1}(s) (1 − p_jk) + P_{j−1}(s − 1) p_jk .
```

The recurrence avoids cancellation between its non-negative summands. This is a stability advantage, not a guarantee
of relative accuracy for arbitrarily small tails: rounding of the inputs, the relative error of 1 − p when p is close
to 1, and the deliberate truncation below come first. (Chen & Liu 1997 survey the Poisson-binomial distribution; the
same recursion is used for credit portfolio losses, Andersen, Sidenius & Basu 2003; electoral-vote totals with margins
of error have been estimated from state-level polls by Kaplan & Barnett 2003.)

**Regions and convolution.** Given (z, mu), the seat distribution of each region is the integral over that region's
factors w_r of its conditional Poisson-binomial distribution; region distributions are then convolved. The national
distribution is the integral of the convolution over (z, mu).

**Events.** A threshold event such as S_k ≥ 326 follows directly from a party marginal. 326 is the majority of the
650-seat House of Commons; since the model covers 631 Great Britain seats, a category with at least 326 of them holds
a majority of the House whatever happens in Northern Ireland. Threshold events of different categories at 326 are
mutually exclusive, so P(no category ≥ 326) = 1 − Σ_k P(S_k ≥ 326) exactly. Because Oth combines different
independents and minor parties, this event means that no model category reaches 326 seats; we report it under that
name (v0.1 called it "hung parliament"). It is not a forecast of government formation or of a working majority.

**Group totals.** A group total (for example Lab + LD) can be computed by summing the group's seat-level win
probabilities before applying the recursion; convolving party marginals would be wrong because the parties compete for
the same seats. A group-total interface is not implemented in this release. Events that constrain two parties
separately (for example Lab ≥ 300 and LD ≥ 60) need a two-dimensional recursion and are not supported.

**Truncation.** Seats whose conditional win probability is within ε = 1e-12 of 0 or 1 are treated as certain. Two
bounds describe the resulting error and should not be confused. The a-priori bound n·ε = 6.3e-10 per party holds for
every value of the factors, because each of the n Bernoulli probabilities moves by at most ε. The values reported per
run (3.1e-12 to 3.3e-12 for the hybrid runs of Table N3) are computed from the probability mass actually dropped,
summed over all eight categories (so each is also a bound for every single party), and are the largest such values over
the nodes evaluated; they say nothing about nodes not evaluated. Neither covers the integration error or
floating-point error.

**Sensitivities.** For a fixed F, ∂P(S_k ≥ m | F)/∂p_ik(F) = P(S_{k,−i} = m − 1 | F), the probability that seat i is
pivotal. Parameter sensitivities follow by integrating P(S_{k,−i} = m − 1 | F) ∂p_ik(F)/∂θ over F (for the mixture,
a weighted sum of the two softmax Jacobians). This identity is stated for completeness; it is not implemented in the
code described here.

**Implementation.** Python with Numba, parallel over factor nodes; exp(eta/τ) is factorised so that the inner loops
only multiply. An active-support variant of the recursion (skipping seats within ε of 0 or 1) gave the same party
PMFs as the dense recursion to 1.4e-14 and was about 2.4 times faster at the same factor points (stored run; see
Appendix F for the re-run script).

## 4. Integration schemes and error components

S2–S4 target the same fixed v0.2.1 model. S1 and the adaptive experiment of Appendix A concern the earlier v0.1
model and are reported separately.

| scheme | national and seat-type factors (9 dims) | regional factors (57 dims) | seat outcomes | used as |
|---|---|---|---|---|
| S1 deterministic tensor (v0.1 engine, v0.1 model) | rotated national axes: trapezoid on the leading axis, Gauss–Hermite on the rest (no seat-type factors in v0.1) | Gauss–Hermite tensor per region | recursion | Appendix A only |
| S2 hybrid (v0.2.1 engine) | scrambled Sobol, 2^m points × R independent scramblings | Gauss–Hermite tensor per region, rule `reg_ns` (points on the regional axes ranked 1–6) | recursion | main engine |
| S3 conditional RQMC | scrambled Sobol over all 66 dimensions | (sampled with the rest) | recursion | reference |
| S4 generative Monte Carlo | sampled | sampled | sampled (Gumbel-max, one scale per seat) | independent check of the implementation |

Write θ for an event probability under the model, θ_H for the value the hybrid scheme converges to with its regional
rule fixed (the outer integral done exactly), and θ̂_H for one hybrid estimate. Then

```
θ̂_H − θ = (θ̂_H − θ_H) + (θ_H − θ) .
```

The scramble standard error estimates the variability of the first term only; it does not address the second, the
mean bias of the regional rule. (This does not mean that the rule error is the same constant at every outer point.)
In more detail, the error of S2 has three components:

1. **Outer randomisation error** (first term). Estimated by the spread over independent scramblings; Owen (1997)
   gives the variance theory of scrambled nets. With few scramblings this estimate can be unstable, in particular when
   rare outer regions contribute (Table N4).
2. **Regional-rule bias** (second term). Deterministic given the rule; it is common to all scramblings and therefore
   invisible in their spread. Changing the rule at the same outer points gives paired differences that are a
   sensitivity measure, not an error bound: Gauss–Hermite rules of different orders are not nested, and the
   differences need not shrink monotonically (Appendix B).
3. **Truncation and floating point.** Bounded as in Section 3 for truncation; negligible here.

A reference from S3 is an estimate with its own sampling error, not the true value. We therefore audit with four
pieces: (a) paired regional-rule ladders at fixed outer points; (b) a conditional RQMC reference; (c) a direct
estimate of the regional-rule bias, comparing the hybrid core with regional-only conditional RQMC at the same outer
points, using independent outer scrambled nets so that the replicate spread is a valid standard error; and (d) an
independent generative Monte Carlo that shares no integration code with the engine.

## 5. Results: two failures of the error bar

Eighteen event probabilities are reported throughout (a Labour majority, seat thresholds for each party, the actual
outcome or more/fewer, no category ≥ 326). They share the same distributions and are not independent: for example a
Labour majority and "no category ≥ 326" are close to complements (correlation −0.9997 across reference replicates), so
the number of events does not add up evidence. Table M1 collects the two failures; the full 18-event tables of the
intermediate steps are in the appendices, while the final comparison (Table N7a) and Table N8 stay here.

{{M1}}

### 5.1 Detection by the generative Monte Carlo

With the default regional rule (3,3,2,2,2,2), two of the 18 events differed from the generative Monte Carlo by about
three standard errors (Table N2 in Appendix D: Grn ≥ 4 z = +3.1, Oth ≥ 5 z = −3.0), while the hybrid's own scramble
standard error for Green was 3.5e-5. The two causes turned out to be different.

### 5.2 Failure 1: the regional rule (Green)

At fixed outer points, raising regional axes 3–6 from 2 to 3 points moved Green by −3.4e-3 and Plaid Cymru by
−5.3e-4; raising axes 4–6 further from 3 to 5 points moved Green again by +6.9e-4 while the other 17 events moved by at
most 4.4e-5 (Table N3, Appendix B). The changes are not monotone. Axis by axis, axis 4 accounts for most of the Green
change (+5.46e-4 from 3 to 9 points, 1e-6 between 7 and 9 points), axis 5 adds +1.6e-4 to +1.7e-4, axis 6 (present
only in Scotland and Wales) does not matter, and the effects of axes 4 and 5 add up (Table N6). Against the
conditional RQMC reference, the default rule is off by +3.2e-3 for Green with a scramble standard error of 3.5e-5
(Table N5 in Appendix D; +3.0e-3 against the pooled reference of Table M1). The outer error bar is therefore not informative about
this error, as the decomposition of Section 4 says it cannot be.

### 5.3 Failure 2: a rare outer contribution (Other ≥ 5)

{{N4}}

Other ≥ 5 hardly moves under the regional rules tried (0.000339 in all of them), so this is not a regional-rule
problem. With four scramblings its estimate was 0.000339 with a standard error of 1.1e-5; larger outer designs gave
0.00041 to 0.00050, consistent with each other and with the reference. Four scramblings had not revealed rare large
outer contributions, so both the estimate and its standard error were off in this case. We do not attach normal-theory
confidence statements to a handful of RQMC replicates and do not generalise this to a claim that standard errors from
some fixed small number of scramblings are always invalid; Jain, Hickernell, Owen & Sorokin (2025, preprint) study
intervals for RQMC with bounded integrands that do not rely on normal theory. We have not implemented them, so no
finite-sample guarantee follows for this audit.

### 5.4 The audited rule

We adopted the rule (3,3,3,9,5,3) and compared it with the pooled references under rules fixed before the additional
runs: three independent conditional RQMC references pooled with equal weight; the hybrid level taken from two
independent outer runs of 2¹² × 8 points; the paired increment from (3,3,3,3,3,3) to the audited rule added; per-event
nominal 95% t approximate intervals with Welch–Satterthwaite degrees of freedom; tolerance 1e-3 (the forecast reports
event probabilities to 0.1 percentage point), with 5e-4 shown for reference. The tolerance was chosen after point
estimates had been seen but before valid bias standard errors and the additional runs existed. This is an exploratory
numerical audit, not a pre-registered experiment.

{{N7a}}

For all 18 events the per-event nominal 95% interval for the difference lies within [−1e-3, 1e-3], and for 14 within
[−5e-4, 5e-4]. This is an exploratory, per-event comparison, not a simultaneous or deterministic error guarantee.
Three intervals exclude zero (Con ≥ 121, Con ≤ 50, Grn ≥ 4), all below 3.5e-4, so small residual biases of that order
are possible. For Green, about +1.1e-4 of the +2.35e-4 difference is the estimated regional-rule bias (Appendix C).
Subtracting the bias estimate gives a secondary table in which 16 of 18 intervals lie within [−1e-3, 1e-3]; the two
others (LD ≥ 50, Ref ≥ 50) are not violations but intervals too wide to confirm the tolerance, because the bias
estimate adds variance.

The tolerance is absolute. For the two rarest events it is larger than the probability itself: Lab ≤ 250 (reference
0.000460) and Oth ≥ 5 (0.000466). Both differences are about +4.3e-5, or 9% of the reference probability, and the upper
ends of their intervals, 1.06e-4 and 9.8e-5, are 23% and 21% of it. For these two events the 18/18 result therefore
does not establish high relative accuracy: relative errors of about 20% are compatible with the comparison (up to 28%
for Oth ≥ 5 in the direct run of Table N7d below). For Con ≥ 250 (0.0045) the tolerance is 22% of the probability and
the interval ends lie within 2.6% of it.

**Direct check at fresh outer points (v0.2.2).** The audited rule was chosen with information from the seed-0 outer
points, which also supply the paired increment in Table N7a. After review, the rule was held fixed and run directly at
new, independent outer points (2¹² × 8, seed 3000, not used before), without the increment, and compared with the same
pooled reference. The design was written down before the run, with no pass/fail rule attached (Table N7d). All 18
intervals again lie within [−1e-3, 1e-3], and 15 within [−5e-4, 5e-4]. Two exclude zero. For Grn ≥ 4 the difference
is +2.4e-4, as in Table N7a. For Con ≥ 121 it is +4.9e-4 (+3.3e-4 in N7a), and the regional-bias estimate for this
event (−4e-5 ± 2.0e-4, Appendix C) does not account for it. A residual difference of a few 1e-4 between the hybrid and
the reference therefore cannot be excluded for Con ≥ 121; it lies inside the 1e-3 tolerance but not inside 5e-4. For
Con ≥ 250 the direct run (0.004497) and the hybrid estimate of N7a (0.004597) differ by about three combined scramble
standard errors, and the reference (0.004545) lies between them. With the seed-0 run (0.004478), the outer runs spread
more than their scramble standard errors suggest, a milder case of the behaviour in Section 5.3.

{{N7d}}

The code's default regional rule remains the fast (3,3,2,2,2,2) so that stored runs reproduce; results in this note
pass the audited rule explicitly.

### 5.5 Conditional RQMC versus conditional iid sampling

{{N8}}

At this design and sample size, the observed iid-to-RQMC variance ratios ranged from 2.9 to 34 for 15 events and from
0.9 to 1.1 for the remaining three (Lab ≤ 250, Ref ≥ 50, Other ≥ 5). These exploratory ratios are based on 16
replicates per method; we do not attach significance statements to them. Ref ≥ 50 has probability about 0.33, so the
absence of a gain is not only a rare-event effect; effective dimension and the placement of factors in the Sobol
sequence may also matter. Owen (1997) gives the asymptotic variance theory of scrambled nets; it is not a guarantee
of a variance reduction for every event at this sample size.

## 6. Practical audit procedure

1. Report the integration scheme for each error component separately. For a hybrid scheme, the scramble standard
   error covers only the outer randomisation.
2. Audit the inner rule at fixed outer points (paired differences) and against a reference that does not use the inner
   rule. Treat rule differences as sensitivities, not bounds.
3. Use enough independent scramblings, and check rare events with larger outer designs; a small number of scramblings
   can understate both the estimate and its error.
4. When the bias of the inner rule is estimated by replicates, use independent outer nets and separate random streams
   for the outer and inner draws (Appendix C).
5. Keep an independent generative simulator that shares no integration code with the engine. In this project it
   flagged each of the allocation errors described here before a reference was available.
6. Fix pooling rules and tolerances before adding runs, and say when they were fixed.

## 7. Limitations

- One model, one election, one set of inputs. The 18 events are dependent.
- The references are themselves estimates. The final comparison is exploratory, with per-event nominal intervals and a
  tolerance chosen after seeing point estimates. The tolerance is absolute; for the two rarest events it exceeds the
  probability itself (Section 5.4).
- Timings are from a shared two-core machine; no performance claims are made. Importance sampling (Glasserman & Li
  2005 for the credit analogue) was not benchmarked.
- Only party marginals are implemented. Group totals follow from the same recursion (Section 3) but have no interface
  in this release; joint constraints on two parties are not supported.
- The model's substantive shortcomings (seat efficiency in 2024, the meaning of the seat-type prior) are outside this
  note and are studied separately.

## 8. Reproducibility

Code, data and stored outputs accompany this note (repository root); the root `README.md` separates three levels: (i)
regenerating this document and its tables from stored outputs, (ii) re-aggregating the stored audit
(`scripts/r1_note_tables_v02.py`), and (iii) re-running the experiments, which takes hours on two cores. Every table
is generated from stored outputs by `scripts/seatdist_note_tables_v0_2.py`, and this document is assembled by
`scripts/build_seatdist_note_v0_2_2.py` (Table N7d by `scripts/seatdist_note_tables_v0_2_2.py`). The random-stream
history is in Appendix C, and run lineage and re-running caveats are in Appendix F; package versions are in
`requirements.txt`. Code is released under the MIT licence; this note and the project's own stored outputs under CC BY
4.0 (attribution: Nari Han); bundled data keep their original licences (`DATA_LICENSES.md`).

## 9. AI usage disclosure

Generative AI tools were used extensively. Claude (Anthropic; Claude Opus 5.5, October 2026) wrote most of the code,
ran the numerical experiments and drafted this note, working from the author's instructions. GPT 6.1 Sol (OpenAI)
reviewed the computations and drafts in several rounds and independently re-ran parts of the audit; changes made in
response are part of the audit record. The author chose the questions, directed and reviewed the work, decided which
results to report and how to word them, and takes full responsibility for the content.

## References

- Andersen, L., Sidenius, J. and Basu, S. (2003). All your hedges in one basket. *Risk*, November 2003, 67–72.
  https://www.risk.net/derivatives/credit-derivatives/1530334/all-your-hedges-one-basket
- British Polling Council (2024). The performance of the polls in the 2024 general election.
  https://www.britishpollingcouncil.org/the-performance-of-the-polls-in-the-2024-general-election/
- Chen, S. X. and Liu, J. S. (1997). Statistical applications of the Poisson-binomial and conditional Bernoulli
  distributions. *Statistica Sinica* 7(4), 875–892. https://www3.stat.sinica.edu.tw/statistica/j7n4/j7n44/j7n44.htm
- Gerstner, T. and Griebel, M. (2003). Dimension-adaptive tensor-product quadrature. *Computing* 71(1), 65–87.
  doi:10.1007/s00607-003-0015-5.
- Glasserman, P. and Li, J. (2005). Importance sampling for portfolio credit risk. *Management Science* 51(11),
  1643–1656. doi:10.1287/mnsc.1050.0415.
- Hanretty, C. (2017). Areal interpolation and the UK's referendum on EU membership. *Journal of Elections, Public
  Opinion and Parties* 27(4), 466–483. doi:10.1080/17457289.2017.1287081.
- Jain, A., Hickernell, F. J., Owen, A. B. and Sorokin, A. G. (2025). Empirical Bernstein and betting confidence
  intervals for randomized quasi-Monte Carlo. arXiv:2504.18677 (v2, 2026). Preprint.
- Kaplan, E. H. and Barnett, A. (2003). A new approach to estimating the probability of winning the presidency.
  *Operations Research* 51(1), 32–40.
- Owen, A. B. (1997). Monte Carlo variance of scrambled net quadrature. *SIAM Journal on Numerical Analysis* 34(5),
  1884–1910. doi:10.1137/S0036142994277468.
- Owen, A. B. (2022). On dropping the first Sobol' point. In A. Keller (ed.), *Monte Carlo and Quasi-Monte Carlo
  Methods (MCQMC 2020)*, Springer Proceedings in Mathematics & Statistics 387, 71–86. doi:10.1007/978-3-030-98319-2_4.
- Perrett, S., Scott, R., Griffiths, J., Prosser, C., Bailey, J., Evans, G., Fieldhouse, E., Green, J. and Mellon, J.
  (2025). 2024 BES Constituency Results with Census and Candidate Data, version 1.0. University of Manchester.
  Dataset. doi:10.48420/28430672. (v0.1 of this note gave a wrong DOI for this file.)
- Trefethen, L. N. and Weideman, J. A. C. (2014). The exponentially convergent trapezoidal rule. *SIAM Review* 56(3),
  385–458.
- UK Parliament. psephology-datasette: election results data. https://github.com/ukparliament/psephology-datasette

## Appendix A. Model v0.1: deterministic tensor quadrature and the adaptive search (C5)

This appendix concerns the earlier model v0.1 (no seat-type factors) and the S1 engine. Its probabilities are not
comparable with those of the v0.2.1 model.

{{N1}}

The reference grid agrees with 3 × 10⁶ generative draws for 17 of 18 events (|z| ≤ 2.5). The exception, Green ≥ 4
(difference 2.1e-3), came from too few points on the regional axes that carry Green's convex response. A larger earlier
grid (24,705 national nodes, 769 s) had ordered the national axes by raw seat influence; the axis that moves the SNP
ranked fourth and received three points, giving errors of 1.4e-2 (SNP ≤ 9) and 1.6e-2 (Other ≥ 5). A grid ladder that
refined only some axes estimated its own error at 6e-4 and missed these; the independent Monte Carlo found them.
Ordering the axes by influence normalised per party fixed the allocation.

**Target-event dimension-adaptive tensor quadrature (exploratory).** We tried to choose the v0.1 grid automatically.
The variables are the point counts on 11 allocation axes (5 national axes and the regional axes ranked 1–6). From a
coarse start, the search tries one more level on each axis and adopts the axis whose refinement changes the 18 event
probabilities most per second of run time; it stops when the largest change falls below 4e-4 or a soft time budget is
exhausted. The search uses no Monte Carlo; the Monte Carlo column is added afterwards. This is a greedy variant of the
dimension-adaptive idea of Gerstner & Griebel (2003) on full tensor grids, with the target defined by the reported
events.

{{N9}}

The final grid (19,250 nodes, 161 s per run, largest difference from Monte Carlo 1.2e-3) has a smaller maximum
difference than our hand design (22,575 nodes, 479 s, 2.1e-3), but it is not better event by event: it is worse on 5
of the 18 events, and three events differ from the 3 × 10⁶ Monte Carlo by more than 2.5 standard errors (Con ≥ 121
z = −4.1, Ref ≥ 50 z = +4.6, Grn ≥ 4 z = −3.4) against one for the hand design (Grn ≥ 4 z = +7.2). The trapezoid
spacing on the leading national axis stayed at the coarsest value (0.7) throughout, consistent in this case with the
fast convergence of the trapezoid rule for smooth integrands with Gaussian decay (Trefethen & Weideman 2014); we do not
claim that result guarantees the accuracy here. Three weaknesses remain: the search took 6,206 s, about 39 final runs;
between 2,970 and 8,250 nodes the one-step criterion fell while the error against Monte Carlo stayed near 3e-3, which
suggests interactions between axes that a one-axis look-ahead cannot see; and the stopping criterion was never
reached. Whether reusing the grid across many scenarios recovers the search cost is untested.

## Appendix B. Regional-rule ladders at fixed outer points (v0.2.1)

{{N3}}

Raising the first two regional axes from 3 to 5 points moves every event by at most 1.1e-4. The regional rule is the
source of the Green discrepancy, the changes are not monotone, and the rule (3,3,3,3,3,3) is not yet converged for
Green.

{{N6}}

Plaid Cymru is insensitive to axes 4–6. Its difference in Table N5 came from the single outer run used there (scramble
SE 1.3e-4); it is absent from the audited comparison (Table N7a: −8.3e-5, interval containing zero; Table N7b:
regional-rule bias −8.5e-7 ± 2.5e-5).

## Appendix C. Regional-rule bias of the audited rule and the secondary comparison

{{N7b}}

{{N7c}}

**Random-stream history.** The stored bias run seeded the outer net with `default_rng([S, r])` and the inner
scramblings with `default_rng([S, r, p, k])`; numpy's SeedSequence zero-pads short entropy, so the first outer point's
first inner stream coincided with the outer stream. Tables N7b and N7c therefore drop the first point of every net
(seven points per net remain). This is an exceptional treatment forced by the stream collision, not a recommended use
of Sobol' points: dropping the first point can damage the digital-net structure (Owen 2022). It does not by itself
imply that the estimate is biased: each scrambled point is marginally uniform, so the mean of any subset of points is
unbiased, and the independence of the 16 nets keeps the replicate standard error valid; but each replicate is no longer
a full eight-point net. The code now separates the streams.

## Appendix D. Full tables for the default hybrid rule (v0.2.1)

{{N2}}

{{N5}}

Differences between two poll-bias scenarios (λ = 0.5 versus 0) are smooth and precise in both schemes when common
random numbers are used, so we no longer claim that smooth scenario differences are specific to deterministic
integration; the Green difference still inherits the regional-rule bias (−0.1675 against −0.1645).

## Appendix E. Post-review follow-up (not part of the reviewed audit)

After the review of v0.1, a follow-up of C5 for the v0.2.1 engine (automatic choice of the regional rule at fixed
outer points, with per-event tolerances, also compared at poll-bias values not used in the search) and a check of
scramble-count choices were run. The follow-up failed its pre-stated success criteria; its later observations are
exploratory. It is reported in a separate C5 follow-up report, which is not part of this bundle, and it does not change
any number or claim in this note.

## Appendix F. Run lineage and re-running caveats

| table | runs (scripts in `scripts/` → outputs in `out/`) |
|---|---|
| M1 | compiled from the outputs of N2, N4 and N7a |
| N1 | `ref_exact2_v0_1.py`, `ref_exact_v0_1.py`, `mc_compare_v0_1.py` → `ref_exact2_v0_1.json`, `ref_exact_v0_1.json`, `mc_compare_std_v0_1.json` |
| N2 | `v02_mc_check.py` → `v02_mc_check.json` |
| N3, N4 | `reg_convergence_v02.py`, `reg_convergence_v02b.py` → `reg_convergence_v02.json`, `reg_convergence_v02b.json` |
| N5 | `cond_mc_v02.py` (reference, seed 20261022), `hybrid_lam_crn_v02.py`; also reads `reg_convergence_v02.json` → `cond_rqmc_N8192_R32_s20261022.json`, `hybrid_lam_crn_v02.json` |
| N6 | `r1_followup_v02.py`, then `r1_green_axis4_v02.py`, then `r1_green_axis5_v02.py` (the last two add rows to the same file; running the first alone resets it) → `r1_followup_v02.json`; kernel benchmark of Section 3 alone: `pb_kernel_bench_v0_1.py` → `pb_kernel_bench_v0_1.json` |
| N7 | references: `cond_mc_v02.py` with seed 20261022 (run directly), 20261103 (called by `r1_green_resolve_v02.py`) and 20261104 (called by `r1_extra_runs_v02.py`); hybrid level: `r1_green_resolve_v02.py`, `r1_extra_runs_v02.py`; increment: `r1_followup_v02.json` (N6 chain); regional bias: `r1_fixed_outer_reps_v02.py`; tables: `r1_note_tables_v02.py` → `r1_note_tables_v02.json`; replicate values of the seed-1 level run and of the increment: `r1_raw_reps_v0_1.py` → `r1_raw_reps_v0_1.json` |
| N7d | `audited_fresh_seed_v0_1.py` → `audited_fresh_seed_v0_1.json` (all replicate values); table: `seatdist_note_tables_v0_2_2.py` |
| N8 | `cond_mc_v02.py` → `cond_rqmc_N2048_R16_s20261020.json`, `cond_iid_N2048_R16_s20261021.json` |
| N9 | `adaptive_quad_v0_1.py` → `adaptive_quad_v0_1.json` |

Caveats for re-running. (i) The stored regional-bias run `out/r1_fixed_outer_reps_v02.json` was produced before the
random-stream fix (Appendix C); the corrected script writes `out/r1_fixed_outer_reps_v02_seedfix.json` and does not
regenerate the stored file, whose re-aggregation (points 2–8 of each net) is what Tables N7 use. (ii) The stored
conditional runs used the dense recursion kernel (`UKSEAT_PB=dense`); the current default is the active-support
kernel, which agrees to 1.4e-14 in the party PMFs. (iii) Up to v0.2.1, the dense side of the kernel benchmark in
`r1_followup_v02.py` was taken from `pb_all`, which follows `UKSEAT_PB` and has defaulted to the active kernel since
v0.2.1, so a default re-run compared the active kernel with itself (found in review, 2026-10-21). The stored file was
produced while the default was still dense; its PMF difference of 1.4e-14 is between two different kernels. The script
now names `pb_dense` explicitly and records each kernel's function name; re-running it also re-runs the N6 audit and
overwrites the stored file, so `pb_kernel_bench_v0_1.py` re-runs the benchmark alone. On the machine used for this
note it gave a maximum PMF difference of 1.37e-14, the stored value, and speed ratios of 2.3 (one λ) and 2.6 (two λ),
against 2.4 and 2.4 stored. (iv) The replicate values behind two inputs of Table N7 were not stored in v0.2.1 (only
means and standard errors). `r1_raw_reps_v0_1.py` re-runs both with their original seeds: the seed-1 level run (2¹² ×
8, rule (3,3,3,3,3,3)) and the audited rule at the seed-0 points (2¹¹ × 4), paired with the stored (3,3,3,3,3,3)
replicates. The means and standard errors reproduce the stored values exactly, and all replicate values and paired
differences are now stored in `r1_raw_reps_v0_1.json`.

Code fix in v0.2.2. In v0.2.1, `QEngine2` and `FastEngine` inherited `ExactEngine.seat_marginals()` (seat-level win
probabilities), which omits the seat-type term L·μ and, for these classes, evaluates the national factor at one node
only; it raised no error. For `QEngine2` with the national factor off, μ fixed at its prior mean and one regional node,
the sum of the returned probabilities differed from the expected seats of `run()` by up to 4.8 seats (Reform 23.08
against 27.87), and by up to 33 seats (Conservatives) with μ at its realised 2024 values; for `FastEngine` with a
five-dimensional outer rule the difference was 4.3 seats. No table or number in this note used the method:
the scripts that call it use `ExactEngine`, for which it is correct. v0.2.2 overrides it in both classes so that it
uses the same outer points (or outer rule), seat-type factors and regional rule as `run()`.

Unit tests: `python tests/test_engine_v0_1.py` (enumeration of small elections against the recursion, truncation
bound, quadrature moments, closed-form mixture against Gumbel-max sampling) and `python tests/test_engine2_v0_2.py`
(the v0.2 core with the seat-type term against enumeration; for `QEngine2` and `FastEngine`, the sum of
`seat_marginals()` against the expected seats of `run()` with the same arguments, scrambling by scrambling, for
national factors on and off, μ fixed or drawn from a normal or Student-t prior, and the default and audited regional
rules; largest difference 1.8e-12 seats).

Data: UK Parliament / House of Commons Library election results (psephology datasette, notional 2019 results on 2024
boundaries and 2024 results; Open Parliament Licence); the BES 2024 constituency results file (Perrett et al. 2025,
doi:10.48420/28430672; CC BY 4.0 according to the University of Manchester record), from which the 2021 census shares
and Chris Hanretty's EU referendum estimates (Hanretty 2017) on 2024 boundaries are taken; British Polling Council
summary of final 2024 polls. v0.1 of this note gave the BES file a DOI that belongs to another dataset. Licences of all
bundled files are listed in `DATA_LICENSES.md`.
