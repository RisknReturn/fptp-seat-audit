### Table M1. The two failures of the default hybrid error bar, and the audited result

Default hybrid: regional rule (3,3,2,2,2,2), 2¹¹ × 4 scrambled Sobol (seed 0); generative MC: 2 × 10⁵ draws (binomial SE); pooled reference and audited-rule difference as in Table N7a (per-event nominal 95% interval).

| event | default hybrid ± scramble SE | generative MC ± SE | pooled reference ± SE | default − reference | audited rule − reference [95% interval] |
|---|---:|---:|---:|---:|---:|
| Grn ≥ 4 | 0.471976 ± 3.5e-5 | 0.46850 ± 1.1e-3 | 0.468985 ± 7.3e-5 | +3.0e-3 | +2.4e-4 [+8.1e-5, +3.9e-4] |
| Oth ≥ 5 | 0.000339 ± 1.1e-5 | 0.00049 ± 4.9e-5 | 0.000466 ± 7.3e-6 | −1.3e-4 | +4.3e-5 [−1.3e-5, +9.8e-5] |

### Table N1. Model v0.1: deterministic grid versus generative Monte Carlo (3 × 10⁶ draws)

Reference grid: national 22,575 nodes (trapezoid h = 0.35 on the leading rotated axis, Gauss–Hermite 7/5/5/3), regional GH (5,5,3,2,2,2); 479 s on 2 shared cores. Earlier reference grid: national axes ordered by raw seat influence (24,705 national nodes, 769 s).

| event | deterministic grid | MC 3 × 10⁶ | z | earlier grid | z |
|---|---:|---:|---:|---:|---:|
| Lab ≥ 326 (majority) | 0.95510 | 0.95493 | +1.3 | 0.95502 | +0.7 |
| Lab ≥ 400 | 0.48515 | 0.48469 | +1.6 | 0.48444 | −0.9 |
| Lab ≥ 450 | 0.11063 | 0.11050 | +0.7 | 0.11047 | −0.2 |
| Lab ≤ 411 | 0.62455 | 0.62506 | −1.8 | 0.62521 | +0.5 |
| Lab ≤ 250 | 0.00028 | 0.00030 | −1.3 | 0.00029 | −1.1 |
| Con ≥ 121 | 0.48408 | 0.48448 | −1.4 | 0.48364 | −2.9 |
| Con < 100 | 0.34571 | 0.34562 | +0.4 | 0.34601 | +1.4 |
| Con ≤ 50 | 0.05496 | 0.05483 | +0.9 | 0.05495 | +0.9 |
| Con ≥ 250 | 0.00397 | 0.00394 | +0.8 | 0.00395 | +0.4 |
| LD ≥ 72 | 0.11925 | 0.11890 | +1.9 | 0.11878 | −0.6 |
| LD ≥ 50 | 0.43848 | 0.43777 | +2.5 | 0.43837 | +2.1 |
| Ref ≥ 5 | 0.96688 | 0.96690 | −0.2 | 0.96676 | −1.3 |
| Ref ≥ 50 | 0.19635 | 0.19612 | +1.0 | 0.19720 | +4.7 |
| Grn ≥ 4 | 0.47418 | 0.47212 | +7.2 | 0.46889 | −11.2 |
| SNP ≤ 9 | 0.09627 | 0.09635 | −0.5 | 0.11036 | +82.3 |
| PC ≥ 4 | 0.45288 | 0.45230 | +2.0 | 0.45252 | +0.8 |
| Oth ≥ 5 | 0.34334 | 0.34374 | −1.4 | 0.32784 | −58.0 |
| No category ≥ 326 | 0.04490 | 0.04506 | −1.3 | 0.04498 | −0.7 |

### Table N2. Model v0.2.1, hybrid scheme with the default regional rule (3,3,2,2,2,2) versus generative MC (2 × 10⁵ draws)

| event | hybrid | scramble SE | MC | z |
|---|---:|---:|---:|---:|
| Lab ≥ 326 (majority) | 0.9380 | 4.6e-4 | 0.9374 | +0.9 |
| Lab ≥ 400 | 0.4356 | 1.6e-4 | 0.4356 | +0.0 |
| Lab ≥ 450 | 0.0957 | 1.9e-4 | 0.0957 | +0.1 |
| Lab ≤ 411 | 0.6679 | 3.4e-4 | 0.6669 | +0.8 |
| Lab ≤ 250 | 0.0004 | 4.6e-5 | 0.0005 | −0.6 |
| Con ≥ 121 | 0.4742 | 2.6e-4 | 0.4739 | +0.3 |
| Con < 100 | 0.3582 | 1.4e-4 | 0.3581 | +0.1 |
| Con ≤ 50 | 0.0613 | 6.0e-5 | 0.0612 | +0.1 |
| Con ≥ 250 | 0.0045 | 7.5e-5 | 0.0045 | +0.1 |
| LD ≥ 72 | 0.1173 | 1.6e-4 | 0.1175 | −0.2 |
| LD ≥ 50 | 0.4209 | 1.6e-4 | 0.4212 | −0.2 |
| Ref ≥ 5 | 0.9862 | 4.5e-5 | 0.9863 | −0.4 |
| Ref ≥ 50 | 0.3262 | 4.6e-4 | 0.3242 | +1.8 |
| Grn ≥ 4 | 0.4720 | 3.5e-5 | 0.4685 | +3.1 |
| SNP ≤ 9 | 0.0877 | 4.5e-4 | 0.0872 | +0.7 |
| PC ≥ 4 | 0.4688 | 1.3e-4 | 0.4676 | +1.1 |
| Oth ≥ 5 | 0.0003 | 1.1e-5 | 0.0005 | −3.0 |
| No category ≥ 326 | 0.0620 | 4.6e-4 | 0.0626 | −0.9 |

### Table N3. Regional rules at fixed outer points (2¹¹ × 4 scrambled Sobol, seed 0): paired differences

Rule = Gauss–Hermite points on the regional factor axes ranked 1–6. All runs share the outer points, so differences isolate the regional rule.

| event | default (3,3,2,2,2,2) | (3,3,3,3,3,3) − default | (5,5,2,2,2,2) − default | (5,5,3,3,3,3) − default | (5,5,5,3,3,3) − default | (3,3,3,5,5,5) − (3,3,3,3,3,3) [paired SE] | outer scramble SE (default) |
|---|---:|---:|---:|---:|---:|---:|---:|
| Lab ≥ 326 (majority) | 0.938037 | −1.5e-4 | −9.0e-7 | −1.6e-4 | −1.6e-4 | +5.9e-6 [6.8e-7] | 4.6e-4 |
| Lab ≥ 400 | 0.435616 | −3.1e-4 | −2.6e-6 | −3.1e-4 | −3.1e-4 | −3.8e-5 [3.2e-7] | 1.6e-4 |
| Lab ≥ 450 | 0.095724 | −1.2e-4 | −5.5e-6 | −1.3e-4 | −1.3e-4 | −2.8e-5 [1.8e-7] | 1.9e-4 |
| Lab ≤ 411 | 0.667883 | +2.8e-4 | +7.4e-6 | +2.8e-4 | +2.8e-4 | +4.0e-5 [6.3e-7] | 3.4e-4 |
| Lab ≤ 250 | 0.000418 | +2.3e-6 | −2.3e-7 | +2.0e-6 | +2.1e-6 | −2.5e-7 [7.9e-8] | 4.6e-5 |
| Con ≥ 121 | 0.474205 | −2.1e-4 | −4.7e-5 | −2.5e-4 | −2.5e-4 | −4.4e-5 [5.2e-7] | 2.6e-4 |
| Con < 100 | 0.358198 | +1.7e-4 | +1.1e-4 | +2.8e-4 | +2.8e-4 | +4.4e-5 [2.5e-7] | 1.4e-4 |
| Con ≤ 50 | 0.061289 | +3.2e-5 | +1.1e-4 | +1.4e-4 | +1.4e-4 | +1.7e-5 [4.2e-7] | 6.0e-5 |
| Con ≥ 250 | 0.004487 | −8.1e-6 | +5.6e-6 | −2.5e-6 | −2.5e-6 | −1.0e-7 [2.8e-8] | 7.5e-5 |
| LD ≥ 72 | 0.117346 | −7.3e-7 | +1.5e-5 | +1.1e-5 | +4.7e-6 | −4.8e-8 [1.0e-7] | 1.6e-4 |
| LD ≥ 50 | 0.420926 | −3.1e-5 | −3.8e-5 | −6.5e-5 | −5.7e-5 | −8.7e-7 [7.0e-8] | 1.6e-4 |
| Ref ≥ 5 | 0.986240 | +2.6e-5 | +9.2e-5 | +1.2e-4 | +1.2e-4 | +3.8e-7 [1.0e-7] | 4.5e-5 |
| Ref ≥ 50 | 0.326199 | +3.3e-4 | −7.6e-5 | +2.5e-4 | +2.5e-4 | −1.1e-5 [4.0e-7] | 4.6e-4 |
| Grn ≥ 4 | 0.471976 | −3.4e-3 | −6.0e-5 | −3.5e-3 | −3.5e-3 | +6.9e-4 [6.1e-6] | 3.5e-5 |
| SNP ≤ 9 | 0.087736 | +1.9e-4 | −5.4e-6 | +1.9e-4 | +1.9e-4 | +4.0e-6 [8.5e-7] | 4.5e-4 |
| PC ≥ 4 | 0.468849 | −5.3e-4 | −2.6e-6 | −5.4e-4 | −5.4e-4 | +9.7e-6 [5.2e-7] | 1.3e-4 |
| Oth ≥ 5 | 0.000339 | −2.8e-7 | +1.8e-8 | −2.6e-7 | −2.6e-7 | −1.9e-8 [6.3e-9] | 1.1e-5 |
| No category ≥ 326 | 0.061953 | +1.5e-4 | +8.5e-7 | +1.6e-4 | +1.6e-4 | −5.9e-6 [6.8e-7] | 4.6e-4 |

Time (2 shared cores): (3,3,2,2,2,2) 47 s, (3,3,3,3,3,3) 133 s, (5,5,2,2,2,2) 100 s, (5,5,3,3,3,3) 362 s, (5,5,5,3,3,3) 585 s. Node-wise Poisson-binomial truncation bound (TV, max over nodes): 3.1e-12 to 3.3e-12.

### Table N4. Outer resolution of a rare event, P(Oth ≥ 5), default regional rule

| outer design | factor points | estimate | scramble SE | replicate range |
|---|---:|---:|---:|---|
| 2¹¹ × 4, seed 0 | 8,192 | 0.000339 | 1.1e-5 | – |
| 2^12 × 8, seed 0 | 32,768 | 0.000410 | 5.1e-5 | 0.000277 – 0.000737 |
| 2^12 × 16, seed 100 | 65,536 | 0.000503 | 2.9e-5 | 0.000366 – 0.000751 |
| 2^13 × 16, seed 200 | 131,072 | 0.000455 | 1.6e-5 | 0.000305 – 0.000539 |

Conditional RQMC reference (no regional quadrature, 8,192 × 32 scrambles): 0.000453 ± 1.1e-5 (estimated SE).

### Table N5. Hybrid rules versus the conditional RQMC reference (8,192 points × 32 scrambles, 66 factor dimensions, all sampled)

Single runs; a difference contains the sampling error of both the hybrid run and the reference. It is not an expected error or an RMSE.

| event | reference ± est. SE | default − ref | default scramble SE | (3,3,3,3,3,3) − ref | λ-difference (0.5 − 0): hybrid ± SE | λ-difference: reference ± SE |
|---|---|---:|---:|---:|---:|---:|
| Lab ≥ 326 (majority) | 0.93784 ± 1.4e-4 | +2.0e-4 | 4.6e-4 | +4.6e-5 | −0.1106 ± 4e-4 | −0.1096 ± 3e-4 |
| Lab ≥ 400 | 0.43551 ± 2.5e-4 | +1.1e-4 | 1.6e-4 | −2.0e-4 | −0.2083 ± 1e-4 | −0.2079 ± 4e-4 |
| Lab ≥ 450 | 0.09554 ± 1.3e-4 | +1.8e-4 | 1.9e-4 | +6.1e-5 | −0.0648 ± 2e-4 | −0.0649 ± 2e-4 |
| Lab ≤ 411 | 0.66717 ± 1.8e-4 | +7.1e-4 | 3.4e-4 | +9.9e-4 | +0.1773 ± 4e-4 | +0.1778 ± 3e-4 |
| Lab ≤ 250 | 0.00044 ± 2.2e-5 | −2.0e-5 | 4.6e-5 | −1.7e-5 | +0.0028 ± 1e-4 | +0.0027 ± 6e-5 |
| Con ≥ 121 | 0.47370 ± 2.1e-4 | +5.0e-4 | 2.6e-4 | +2.9e-4 | +0.2786 ± 5e-4 | +0.2792 ± 4e-4 |
| Con < 100 | 0.35879 ± 2.3e-4 | −5.9e-4 | 1.4e-4 | −4.2e-4 | −0.2242 ± 2e-4 | −0.2246 ± 3e-4 |
| Con ≤ 50 | 0.06151 ± 1.6e-4 | −2.2e-4 | 6.0e-5 | −1.9e-4 | −0.0499 ± 2e-5 | −0.0500 ± 2e-4 |
| Con ≥ 250 | 0.00457 ± 5.7e-5 | −8.0e-5 | 7.5e-5 | −8.8e-5 | +0.0282 ± 1e-4 | +0.0282 ± 2e-4 |
| LD ≥ 72 | 0.11720 ± 2.0e-4 | +1.5e-4 | 1.6e-4 | +1.5e-4 | −0.0612 ± 2e-4 | −0.0610 ± 3e-4 |
| LD ≥ 50 | 0.42141 ± 2.9e-4 | −4.8e-4 | 1.6e-4 | −5.1e-4 | −0.1546 ± 2e-4 | −0.1552 ± 4e-4 |
| Ref ≥ 5 | 0.98613 ± 6.0e-5 | +1.1e-4 | 4.5e-5 | +1.3e-4 | −0.0047 ± 2e-6 | −0.0047 ± 2e-5 |
| Ref ≥ 50 | 0.32614 ± 3.2e-4 | +6.0e-5 | 4.6e-4 | +3.9e-4 | −0.0375 ± 3e-5 | −0.0376 ± 1e-4 |
| Grn ≥ 4 | 0.46881 ± 1.3e-4 | +3.2e-3 | 3.5e-5 | −2.5e-4 | −0.1675 ± 6e-5 | −0.1645 ± 9e-5 |
| SNP ≤ 9 | 0.08800 ± 1.8e-4 | −2.6e-4 | 4.5e-4 | −6.7e-5 | +0.0079 ± 4e-5 | +0.0079 ± 2e-5 |
| PC ≥ 4 | 0.46804 ± 5.9e-5 | +8.1e-4 | 1.3e-4 | +2.8e-4 | +0.0381 ± 5e-5 | +0.0380 ± 2e-5 |
| Oth ≥ 5 | 0.00045 ± 1.1e-5 | −1.1e-4 | 1.1e-5 | −1.1e-4 | +0.0001 ± 8e-6 | +0.0001 ± 2e-6 |
| No category ≥ 326 | 0.06215 ± 1.3e-4 | −2.0e-4 | 4.6e-4 | −4.4e-5 | +0.1103 ± 4e-4 | +0.1094 ± 3e-4 |

Measured times (2 shared cores, indicative only): default hybrid 47 s for one λ in the rule ladder, 69 s and 79 s in the two λ runs of the CRN comparison; (3,3,3,3,3,3) 133 s; reference 1336 s (dense recursion kernel).

### Table N6. Green: axis-by-axis regional audit at the same outer points (baseline (3,3,3,3,3,3))

| rule | Grn ≥ 4 − baseline [paired SE] | PC ≥ 4 − baseline | Grn ≥ 4 | time |
|---|---:|---:|---:|---:|
| (3,3,3,5,3,3) | +5.198e-4 [6.1e-6] | +1.2e-5 | 0.469075 | 227 s |
| (3,3,3,3,5,3) | +1.690e-4 [6.8e-7] | −2.3e-6 | 0.468724 | 180 s |
| (3,3,3,3,3,5) | +1.388e-17 [2.7e-17] | −1.7e-14 | 0.468555 | 130 s |
| (3,3,3,7,3,3) | +5.454e-4 [7.1e-6] | +1.0e-5 | 0.469101 | 249 s |
| (3,3,3,9,3,3) | +5.464e-4 [7.3e-6] | +1.0e-5 | 0.469102 | 318 s |
| (3,3,3,3,7,3) | +1.607e-4 [4.9e-7] | −2.4e-6 | 0.468716 | 240 s |
| (3,3,3,9,5,3) | +7.154e-4 [7.1e-6] | +7.9e-6 | 0.469271 | 516 s |

### Table N7a. Audited rule (3,3,3,9,5,3) minus the pooled references (three conditional RQMC references, 96 scrambles)

Per-event nominal 95% t approximate intervals (Welch–Satterthwaite df); not simultaneous over the 18 events.

| event | pooled reference ± SE | difference | 95% interval (df) | within 1e-3 | within 5e-4 |
|---|---|---:|---:|:-:|:-:|
| Lab ≥ 326 (majority) | 0.937714 ± 9.5e-5 | −3.9e-5 | [−3.4e-4, +2.6e-4] (38) | yes | yes |
| Lab ≥ 400 | 0.435419 ± 1.3e-4 | −2.3e-4 | [−5.6e-4, +1.0e-4] (76) | yes | no |
| Lab ≥ 450 | 0.095696 ± 8.2e-5 | +2.4e-4 | [−2.2e-5, +5.0e-4] (37) | yes | no |
| Lab ≤ 411 | 0.667448 ± 1.2e-4 | +3.3e-4 | [−6.5e-5, +7.2e-4] (37) | yes | no |
| Lab ≤ 250 | 0.000460 ± 1.4e-5 | +4.3e-5 | [−2.1e-5, +1.1e-4] (23) | yes | yes |
| Con ≥ 121 | 0.473584 ± 1.2e-4 | +3.3e-4 | [+6.9e-5, +6.0e-4] (103) | yes | no |
| Con < 100 | 0.358526 ± 1.3e-4 | −1.4e-6 | [−2.7e-4, +2.7e-4] (105) | yes | yes |
| Con ≤ 50 | 0.061437 ± 8.9e-5 | −1.9e-4 | [−3.8e-4, −4.6e-6] (109) | yes | yes |
| Con ≥ 250 | 0.004545 ± 2.9e-5 | +5.2e-5 | [−1.6e-5, +1.2e-4] (93) | yes | yes |
| LD ≥ 72 | 0.117263 ± 1.2e-4 | +2.2e-6 | [−2.5e-4, +2.6e-4] (107) | yes | yes |
| LD ≥ 50 | 0.421290 ± 1.6e-4 | −1.0e-4 | [−4.5e-4, +2.4e-4] (109) | yes | yes |
| Ref ≥ 5 | 0.986200 ± 3.8e-5 | −1.8e-5 | [−1.4e-4, +1.0e-4] (39) | yes | yes |
| Ref ≥ 50 | 0.326156 ± 1.5e-4 | +9.9e-5 | [−2.8e-4, +4.8e-4] (80) | yes | yes |
| Grn ≥ 4 | 0.468985 ± 7.3e-5 | +2.4e-4 | [+8.1e-5, +3.9e-4] (111) | yes | yes |
| SNP ≤ 9 | 0.087876 ± 9.1e-5 | +1.6e-4 | [−8.6e-5, +4.0e-4] (59) | yes | yes |
| PC ≥ 4 | 0.468094 ± 4.0e-5 | −8.3e-5 | [−1.8e-4, +1.2e-5] (84) | yes | yes |
| Oth ≥ 5 | 0.000466 ± 7.3e-6 | +4.3e-5 | [−1.3e-5, +9.8e-5] (16) | yes | yes |
| No category ≥ 326 | 0.062279 ± 9.5e-5 | +4.4e-5 | [−2.5e-4, +3.4e-4] (38) | yes | yes |

Within 1e-3: 18/18; within 5e-4: 14/18.

### Table N7c. As N7a after subtracting the estimated regional-rule bias (secondary)

Per-event nominal 95% t approximate intervals (Welch–Satterthwaite df); not simultaneous over the 18 events.

| event | pooled reference ± SE | difference | 95% interval (df) | within 1e-3 | within 5e-4 |
|---|---|---:|---:|:-:|:-:|
| Lab ≥ 326 (majority) | 0.937714 ± 9.5e-5 | −7.6e-6 | [−3.1e-4, +3.0e-4] (42) | yes | yes |
| Lab ≥ 400 | 0.435419 ± 1.3e-4 | −2.1e-4 | [−5.8e-4, +1.6e-4] (91) | yes | no |
| Lab ≥ 450 | 0.095696 ± 8.2e-5 | +2.6e-4 | [−2.0e-5, +5.3e-4] (45) | yes | no |
| Lab ≤ 411 | 0.667448 ± 1.2e-4 | +4.0e-4 | [−4.6e-5, +8.5e-4] (51) | yes | no |
| Lab ≤ 250 | 0.000460 ± 1.4e-5 | +4.0e-5 | [−2.4e-5, +1.0e-4] (23) | yes | yes |
| Con ≥ 121 | 0.473584 ± 1.2e-4 | +3.8e-4 | [−1.1e-4, +8.6e-4] (31) | yes | no |
| Con < 100 | 0.358526 ± 1.3e-4 | +4.4e-4 | [−9.9e-5, +9.7e-4] (28) | yes | no |
| Con ≤ 50 | 0.061437 ± 8.9e-5 | −6.0e-5 | [−4.3e-4, +3.1e-4] (27) | yes | yes |
| Con ≥ 250 | 0.004545 ± 2.9e-5 | +6.6e-5 | [−7.2e-6, +1.4e-4] (108) | yes | yes |
| LD ≥ 72 | 0.117263 ± 1.2e-4 | −2.5e-4 | [−6.5e-4, +1.5e-4] (42) | yes | no |
| LD ≥ 50 | 0.421290 ± 1.6e-4 | +2.5e-4 | [−5.7e-4, +1.1e-3] (23) | **no** | no |
| Ref ≥ 5 | 0.986200 ± 3.8e-5 | +7.3e-5 | [−1.2e-4, +2.6e-4] (36) | yes | yes |
| Ref ≥ 50 | 0.326156 ± 1.5e-4 | +3.3e-4 | [−3.9e-4, +1.0e-3] (29) | **no** | no |
| Grn ≥ 4 | 0.468985 ± 7.3e-5 | +1.3e-4 | [−3.7e-5, +2.9e-4] (125) | yes | yes |
| SNP ≤ 9 | 0.087876 ± 9.1e-5 | +1.9e-4 | [−3.1e-4, +6.9e-4] (26) | yes | no |
| PC ≥ 4 | 0.468094 ± 4.0e-5 | −8.2e-5 | [−1.9e-4, +2.5e-5] (97) | yes | yes |
| Oth ≥ 5 | 0.000466 ± 7.3e-6 | +4.3e-5 | [−1.3e-5, +9.8e-5] (16) | yes | yes |
| No category ≥ 326 | 0.062279 ± 9.5e-5 | +1.3e-5 | [−2.9e-4, +3.2e-4] (42) | yes | yes |

Within 1e-3: 16/18; within 5e-4: 9/18.

### Table N7b. Regional-rule bias of the audited rule: hybrid minus regional-only conditional RQMC at the same outer points

16 independent outer scrambled Sobol nets; mean of points 2–8 of each net (the first point is dropped because its inner random stream coincided with the outer one in the stored run); SE from the 16 replicate means, t₁₅.

| event | bias ± SE | 95% interval |
|---|---:|---:|
| Lab ≥ 326 (majority) | −3.1e-5 ± 3.4e-5 | [−1.0e-4, +4.0e-5] |
| Lab ≥ 400 | −1.9e-5 ± 8.0e-5 | [−1.9e-4, +1.5e-4] |
| Lab ≥ 450 | −1.5e-5 ± 4.4e-5 | [−1.1e-4, +7.9e-5] |
| Lab ≤ 411 | −7.3e-5 ± 1.1e-4 | [−3.1e-4, +1.6e-4] |
| Lab ≤ 250 | +2.4e-6 ± 2.4e-6 | [−2.7e-6, +7.4e-6] |
| Con ≥ 121 | −4.1e-5 ± 2.0e-4 | [−4.6e-4, +3.8e-4] |
| Con < 100 | −4.4e-4 ± 2.2e-4 | [−9.1e-4, +3.6e-5] |
| Con ≤ 50 | −1.3e-4 ± 1.6e-4 | [−4.6e-4, +2.0e-4] |
| Con ≥ 250 | −1.4e-5 ± 1.5e-5 | [−4.6e-5, +1.7e-5] |
| LD ≥ 72 | +2.5e-4 ± 1.5e-4 | [−6.5e-5, +5.7e-4] |
| LD ≥ 50 | −3.5e-4 ± 3.6e-4 | [−1.1e-3, +4.0e-4] |
| Ref ≥ 5 | −9.1e-5 ± 7.2e-5 | [−2.4e-4, +6.2e-5] |
| Ref ≥ 50 | −2.3e-4 ± 3.0e-4 | [−8.6e-4, +4.0e-4] |
| Grn ≥ 4 | +1.1e-4 ± 2.6e-5 | [+5.5e-5, +1.6e-4] |
| SNP ≤ 9 | −3.0e-5 ± 2.1e-4 | [−4.8e-4, +4.2e-4] |
| PC ≥ 4 | −8.5e-7 ± 2.5e-5 | [−5.4e-5, +5.2e-5] |
| Oth ≥ 5 | −3.1e-8 ± 6.1e-8 | [−1.6e-7, +9.8e-8] |
| No category ≥ 326 | +3.1e-5 ± 3.4e-5 | [−4.0e-5, +1.0e-4] |

### Table N8. Conditional sampling: observed variance ratio iid / RQMC at the same N (2,048 points, 16 replicates each)

| event | ratio | event | ratio | event | ratio |
|---|---:|---|---:|---|---:|
| Lab ≥ 326 (majority) | 10.5 | Lab ≥ 400 | 16.0 | Lab ≥ 450 | 5.8 |
| Lab ≤ 411 | 7.0 | Lab ≤ 250 | 0.9 | Con ≥ 121 | 13.1 |
| Con < 100 | 10.1 | Con ≤ 50 | 3.9 | Con ≥ 250 | 2.9 |
| LD ≥ 72 | 6.9 | LD ≥ 50 | 9.9 | Ref ≥ 5 | 3.2 |
| Ref ≥ 50 | 1.1 | Grn ≥ 4 | 33.8 | SNP ≤ 9 | 18.0 |
| PC ≥ 4 | 25.5 | Oth ≥ 5 | 1.0 | No category ≥ 326 | 10.4 |

### Table N9. Target-event dimension-adaptive tensor quadrature (model v0.1): search trace

Greedy rule: refine the axis with the largest change of the 18 event probabilities per second. MC is not used for selection or stopping; the MC column is shown afterwards.

| step | grid (national; regional) | nodes | time per run | max \|grid − MC\| | cumulative search time |
|---:|---|---:|---:|---:|---:|
| 0 | trapezoid h = 0.7, GH 2/2/2/2; regional (2,2,2,2,2,2) | 352 | 1.3 s | 0.0180 | 17 s |
| 1 | trapezoid h = 0.7, GH 2/2/3/2; regional (2,2,2,2,2,2) | 528 | 1.6 s | 0.0218 | 42 s |
| 2 | trapezoid h = 0.7, GH 2/3/3/2; regional (2,2,2,2,2,2) | 792 | 2.2 s | 0.0165 | 81 s |
| 3 | trapezoid h = 0.7, GH 3/3/3/2; regional (2,2,2,2,2,2) | 1,188 | 4.2 s | 0.0056 | 139 s |
| 4 | trapezoid h = 0.7, GH 3/3/3/2; regional (2,3,2,2,2,2) | 1,188 | 4.8 s | 0.0055 | 217 s |
| 5 | trapezoid h = 0.7, GH 3/3/3/3; regional (2,3,2,2,2,2) | 1,782 | 7.2 s | 0.0047 | 426 s |
| 6 | trapezoid h = 0.7, GH 3/3/5/3; regional (2,3,2,2,2,2) | 2,970 | 12.5 s | 0.0029 | 832 s |
| 7 | trapezoid h = 0.7, GH 3/3/5/3; regional (2,3,2,3,2,2) | 2,970 | 15.9 s | 0.0033 | 1,146 s |
| 8 | trapezoid h = 0.7, GH 3/3/5/3; regional (3,3,2,3,2,2) | 2,970 | 23.6 s | 0.0032 | 1,793 s |
| 9 | trapezoid h = 0.7, GH 5/3/5/3; regional (3,3,2,3,2,2) | 4,950 | 46.8 s | 0.0032 | 2,509 s |
| 10 | trapezoid h = 0.7, GH 5/5/5/3; regional (3,3,2,3,2,2) | 8,250 | 80.1 s | 0.0032 | 3,801 s |
| 11 | trapezoid h = 0.7, GH 5/5/5/5; regional (3,3,2,3,2,2) | 13,750 | 118.3 s | 0.0016 | 6,206 s |
| final | trapezoid h = 0.7, GH 7/5/5/5; regional (3,3,2,3,2,2) | 19,250 | 161 s | 0.0012 | – |
| hand design (N1) | trapezoid h = 0.35, GH 7/5/5/3; regional (5,5,3,2,2,2) | 22,575 | 479 s | 0.0021 | – |

