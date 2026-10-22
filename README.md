# fptp-seat-audit

[![DOI](https://zenodo.org/badge/DOI/__CONCEPT_DOI__.svg)](https://doi.org/__CONCEPT_DOI__)

**A numerical audit of factor-based FPTP seat distributions: a Great Britain 2024 case study.**
Code, data and stored outputs for the computation note [`docs/SEATDIST_NOTE_v0_2_2.md`](docs/SEATDIST_NOTE_v0_2_2.md)
(draft v0.2.2, not peer reviewed).

> **Scope.** The model here is a fixed computational case (model v0.2.1, final 2024 polls). It is not a forecasting product,
> and it misjudged Labour's 2024 seat efficiency badly. The subject is how to compute, and how to check, the seat
> distribution implied by a given factor model.

## What the note shows

If seat outcomes are independent once a set of common factors is fixed (national polling error, regional swings,
seat-type factors), each party's seat distribution can be computed for every value of the factors with the
Poisson-binomial recursion familiar from credit portfolio loss models, followed by a convolution across regions.
Given the factors this is algebraically exact; the implementation adds floating-point error and a bounded
support-truncation error, and the integral over the 66 factor dimensions is the dominant approximation.

The scheme used most is a hybrid: scrambled Sobol points over nine national and seat-type factors, and a
Gauss–Hermite rule over the regional factors. Its usual error bar, the spread over independent scramblings, missed two
errors of different kinds:

1. **A regional-rule bias.** P(Green ≥ 4 seats) was off by about 3e-3, while its scramble standard error was 3.5e-5.
   The outer error bar cannot see the mean bias of the inner rule.
2. **An under-resolved rare outer contribution.** With four scramblings, both the estimate of P(Other ≥ 5) and its
   standard error were unstable.

Four computation paths separate the error components:
- paired regional-rule changes at fixed outer points
- a conditional RQMC reference without regional quadrature
- a regional-only RQMC bias estimate
- an independent generative Monte Carlo

With a strengthened regional rule, the per-event nominal 95% intervals for the hybrid–reference differences lie
within [−1e-3, 1e-3] for all 18 reported events (16 of 18 in a secondary, bias-corrected comparison). This is an
exploratory per-event comparison, not a simultaneous or deterministic error guarantee.

The note claims no novelty for the recursion and no general speed advantage over Monte Carlo. Its contribution is a
reproducible case in which two failures of a common error bar are each located by a different computation path.

## Install

```
python -m pip install -r requirements.txt
```

Tested versions are pinned in `requirements.txt` (Python 3.13.16 on Linux x86-64). Levels 1 and 2 below need only
numpy and scipy. Level 3, the example and the tests also need pandas and Numba.

On Windows, run Python in UTF-8 mode (`python -X utf8 ...` or `set PYTHONUTF8=1`). Generated Markdown may then differ
from the stored files only in line endings (CRLF), and generated JSON in the last digits of floating-point values
(differences of order 1e-14 were seen in an independent run on Windows). The tables and all judgements are the same.

## Three levels of reproduction

**1. Regenerate the note and its tables from stored outputs (seconds).**

```
python scripts/seatdist_note_tables_v0_2.py      # out/seatdist_note_tables_v0_2.md (Tables M1, N1-N9)
python scripts/seatdist_note_tables_v0_2_2.py    # out/seatdist_note_tables_v0_2_2.md (Table N7d)
python scripts/build_seatdist_note_v0_2_2.py     # docs/SEATDIST_NOTE_v0_2_2.md
```

`scripts/build_seatdist_note_v0_2_1.py` still rebuilds the previous version, `docs/SEATDIST_NOTE_v0_2_1.md`.

**2. Re-aggregate the stored audit (seconds).** This recomputes the final comparison (Tables N7a–c) from the stored
runs.

```
python scripts/r1_note_tables_v02.py             # out/r1_note_tables_v02.json, then rerun level 1
```

**3. Re-run the experiments (hours on two cores).** The scripts and their order are listed table by table in
Appendix F of the note. Three caveats:
- The stored regional-bias run predates a random-stream fix. The corrected script writes a new file and does not
  overwrite the stored one (Appendix C).
- The stored conditional runs used the dense recursion kernel (`UKSEAT_PB=dense`). The current default kernel agrees
  to 1.4e-14.
- `scripts/r1_followup_v02.py` re-runs the kernel benchmark and the Green axis audit together and overwrites its stored
  file. To re-check the kernel benchmark alone, run `scripts/pb_kernel_bench_v0_1.py`.

The continuous-integration workflow (`.github/workflows/tests.yml`) runs the unit tests and checks that levels 1 and 2
reproduce the committed note and tables.

## Example with the audited regional rule

```
python scripts/example_audited_v0_2.py 8 4       # outer 2^8 x 4 scrambled Sobol, regional rule (3,3,3,9,5,3)
```

The engine's constructor default regional rule `(3,3,2,2,2,2)` is kept only so that stored runs reproduce. It has the
Green ≥ 4 bias of about +3e-3 described above. Pass `REG_NS_AUDITED` from `ukseat/engine2.py` explicitly, as the
example does.

The event stored as `Hung (no party>=326)` means that no model category reaches 326 seats, with Oth (independents and
minor parties together) counted as one category. The note reports it under that name. It is not a forecast of
government formation.

Seat-level win probabilities for the same model and outer points: `QEngine2(...).seat_marginals(...)` takes the
arguments of `run(...)`. It is far slower than `run()` for large regional rules. (In v0.2.1 this method was inherited
and wrong for `QEngine2` and `FastEngine`; see Appendix F of the note.)

## Tests

```
python tests/test_engine_v0_1.py
python tests/test_engine2_v0_2.py
```

The first file checks the recursion against enumeration of small elections, the truncation bound, the quadrature
moments, and the closed-form mixture against Gumbel-max sampling. The second checks the v0.2 core (with the seat-type
term) against enumeration, and that the seat-level probabilities of `QEngine2` and `FastEngine` add up to the expected
seats of `run()` with the same arguments.

## Repository layout

| path | content |
|---|---|
| `docs/SEATDIST_NOTE_v0_2_2.md` | the note (current version) and its template. `v0_1`, `v0_2` and `v0_2_1` are earlier versions, kept as history |
| `ukseat/` | model (`model.py`, `model2.py`), engines (`engine.py`, `engine2.py`, Numba core `core.py`), quadrature, 2024 inputs (`config_2024.py`), events (`report.py`) |
| `scripts/` | data preparation, experiments, table and note generation |
| `data/` | inputs (sources and licences in `DATA_LICENSES.md`) |
| `out/` | stored outputs (JSON, logs, generated tables) |
| `tests/` | unit tests |

## Citation

Cite all versions with the concept DOI [__CONCEPT_DOI__](https://doi.org/__CONCEPT_DOI__), which resolves to the
latest version; version 0.2.2 is [__V022_DOI__](https://doi.org/__V022_DOI__), and the DOI of this version (0.2.3) is
listed on its Zenodo record. GitHub's "Cite this repository" button uses [`CITATION.cff`](CITATION.cff).

> Han, N. (2026). *fptp-seat-audit: a numerical audit of factor-based FPTP seat distributions (Great Britain 2024
> case study)*, version 0.2.3. Risk & Return Co., Ltd. Zenodo. https://doi.org/__CONCEPT_DOI__

## Licences

- **Code** (`ukseat/`, `scripts/`, `tests/`): MIT licence (`LICENSE`).
- **The note and the project's own stored outputs** (`docs/`, `out/`): CC BY 4.0 (attribution: Nari Han).
- **Data:** files keep the terms of their sources. UK Parliament election results are under the Open Parliament
  Licence v3.0 (Contains Parliamentary information licensed under the Open Parliament Licence v3.0); the BES 2024
  constituency results file is under CC BY 4.0. Details and attribution are in [`DATA_LICENSES.md`](DATA_LICENSES.md).

## AI usage disclosure

Generative AI tools were used extensively. Claude (Anthropic; Claude Opus 5.5, October 2026) wrote most of the code,
ran the numerical experiments and drafted the note, working from the author's instructions. GPT 6.1 Sol (OpenAI)
reviewed the computations and drafts in several rounds and independently re-ran parts of the audit; changes made in
response are part of the audit record. The author chose the questions, directed and reviewed the work, decided which
results to report and how to word them, and takes full responsibility for the content.

## Contact

Comments and corrections are welcome through GitHub Issues.
