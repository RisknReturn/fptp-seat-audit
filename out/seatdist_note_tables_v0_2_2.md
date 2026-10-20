### Table N7d. Audited rule (3,3,3,9,5,3) run directly at fresh outer points, minus the pooled reference of Table N7a

Hybrid: rule fixed in advance, 2¹² × 8 scrambled Sobol, seed 3000 (not used before), λ = 0, no paired increment; reference and interval method as in Table N7a (per-event nominal 95% t approximate intervals, Welch–Satterthwaite df; not simultaneous). No pass/fail rule was attached to this run.

| event | pooled reference ± SE | direct run ± scramble SE | difference | 95% interval (df) | N7a difference | within 1e-3 | within 5e-4 |
|---|---|---|---:|---:|---:|:-:|:-:|
| Lab ≥ 326 (majority) | 0.937714 ± 9.5e-5 | 0.937809 ± 1.0e-4 | +9.5e-5 | [−1.9e-4, +3.8e-4] (23) | −3.9e-5 | yes | yes |
| Lab ≥ 400 | 0.435419 ± 1.3e-4 | 0.435313 ± 1.8e-4 | −1.1e-4 | [−5.8e-4, +3.7e-4] (16) | −2.3e-4 | yes | no |
| Lab ≥ 450 | 0.095696 ± 8.2e-5 | 0.095656 ± 1.4e-4 | −4.0e-5 | [−3.8e-4, +3.0e-4] (13) | +2.4e-4 | yes | yes |
| Lab ≤ 411 | 0.667448 ± 1.2e-4 | 0.667853 ± 2.0e-4 | +4.1e-4 | [−1.1e-4, +9.2e-4] (13) | +3.3e-4 | yes | no |
| Lab ≤ 250 | 0.000460 ± 1.4e-5 | 0.000454 ± 3.7e-5 | −5.8e-6 | [−9.5e-5, +8.3e-5] (9) | +4.3e-5 | yes | yes |
| Con ≥ 121 | 0.473584 ± 1.2e-4 | 0.474073 ± 6.9e-5 | +4.9e-4 | [+2.1e-4, +7.6e-4] (67) | +3.3e-4 | yes | no |
| Con < 100 | 0.358526 ± 1.3e-4 | 0.358435 ± 6.0e-5 | −9.1e-5 | [−3.8e-4, +2.0e-4] (89) | −1.4e-6 | yes | yes |
| Con ≤ 50 | 0.061437 ± 8.9e-5 | 0.061320 ± 5.2e-5 | −1.2e-4 | [−3.2e-4, +8.8e-5] (66) | −1.9e-4 | yes | yes |
| Con ≥ 250 | 0.004545 ± 2.9e-5 | 0.004497 ± 2.7e-5 | −4.8e-5 | [−1.3e-4, +3.4e-5] (29) | +5.2e-5 | yes | yes |
| LD ≥ 72 | 0.117263 ± 1.2e-4 | 0.117150 ± 5.3e-5 | −1.1e-4 | [−3.7e-4, +1.4e-4] (88) | +2.2e-6 | yes | yes |
| LD ≥ 50 | 0.421290 ± 1.6e-4 | 0.421327 ± 1.1e-4 | +3.7e-5 | [−3.5e-4, +4.3e-4] (55) | −1.0e-4 | yes | yes |
| Ref ≥ 5 | 0.986200 ± 3.8e-5 | 0.986105 ± 5.4e-5 | −9.5e-5 | [−2.4e-4, +4.6e-5] (16) | −1.8e-5 | yes | yes |
| Ref ≥ 50 | 0.326156 ± 1.5e-4 | 0.326099 ± 1.2e-4 | −5.7e-5 | [−4.5e-4, +3.4e-4] (43) | +9.9e-5 | yes | yes |
| Grn ≥ 4 | 0.468985 ± 7.3e-5 | 0.469229 ± 2.7e-5 | +2.4e-4 | [+9.0e-5, +4.0e-4] (98) | +2.4e-4 | yes | yes |
| SNP ≤ 9 | 0.087876 ± 9.1e-5 | 0.087886 ± 1.4e-4 | +9.7e-6 | [−3.5e-4, +3.7e-4] (14) | +1.6e-4 | yes | yes |
| PC ≥ 4 | 0.468094 ± 4.0e-5 | 0.468090 ± 4.2e-5 | −3.9e-6 | [−1.2e-4, +1.2e-4] (23) | −8.3e-5 | yes | yes |
| Oth ≥ 5 | 0.000466 ± 7.3e-6 | 0.000515 ± 3.5e-5 | +4.9e-5 | [−3.6e-5, +1.3e-4] (8) | +4.3e-5 | yes | yes |
| No category ≥ 326 | 0.062279 ± 9.5e-5 | 0.062186 ± 1.0e-4 | −9.3e-5 | [−3.8e-4, +1.9e-4] (24) | +4.4e-5 | yes | yes |

Within 1e-3: 18/18; within 5e-4: 15/18. Run time 43 min on two cores; largest per-node truncation bound 3.7e-12.
