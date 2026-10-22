# Data sources and licences

The MIT licence in `LICENSE` covers the code (`ukseat/`, `scripts/`, `tests/`). The note and the project's own stored
outputs (`docs/`, `out/`) are released under CC BY 4.0 (attribution: Nari Han). Bundled data files keep the terms of
their sources, listed below; the CC BY 4.0 of the note does not change them.

| file | content | source | licence / terms | checked |
|---|---|---|---|---|
| `data/gb2024_seats_v0_1.csv` | GB seats on 2024 boundaries: 2019 notional and 2024 actual votes | UK Parliament `psephology-datasette` (`psephology.db`), compiled by the House of Commons Library from the UK Parliament election results website, https://github.com/ukparliament/psephology-datasette | Open Parliament Licence, as stated in the repository's `metadata.json` (data licence; the repository has no separate code licence file) | 2026-10-21, `metadata.json` |
| `data/old_bound_hist_v0_1.csv` | 2010–2019 actual votes on 2010 boundaries (calibration) | same as above | Open Parliament Licence (as above) | as above |
| `data/gb2024_v02_covariates.csv` | census shares (Muslim, degree, age, ethnicity), EU referendum Leave estimates, 2024 nominations, derived per seat (columns selected and renamed, age bands summed; `scripts/prep_v02.py`) | Census and Leave columns (`leave`, `muslim`, `degree`, `white_british`, `age20_29`, `age65p`): 2024 BES Constituency Results with Census and Candidate Data (Perrett, Scott, Griffiths, Prosser, Bailey, Evans, Fieldhouse, Green and Mellon 2025; made available 19 February 2025; doi:10.48420/28430672). Its census variables come from the 2021 Census (England and Wales) and, for Scottish constituencies, from Scotland's Census 2022, as compiled in the BES file (its codebook notes that the Scottish census was fielded in 2022 and that not all census information is available for Scotland); its Leave estimates come from Chris Hanretty and are included in the BES file with his permission (cite Hanretty 2017). Nomination columns (`stand24_*`, `ind_sitting_mp`, `wpb24`, `oth_returning_share`) and seat identifiers: UK Parliament `psephology.db`, as in the first row | BES columns: CC BY 4.0, as stated in the University of Manchester record, https://research.manchester.ac.uk/en/datasets/2024-bes-constituency-results-with-census-and-candidate-data/ . This applies to this results dataset only, not to the BES panel survey files. Nomination columns: Open Parliament Licence (as above) | 2026-10-21, University of Manchester record |
| `data/bes/` (original BES results file and codebook) | the source file of the row above | as above | CC BY 4.0 (as above). Not included in the bundle because nothing but `prep_v02.py` needs it; download it from the BES or University of Manchester page to rebuild the covariates file. | as above |
| `data/calib_v0_1.json`, `data/calib_mu_v02.json`, `data/realised_mu_2024.json`, `data/regional_dev_v0_1.csv`, `data/seat_resid_v0_1.csv` | parameters and residuals computed by this project's scripts from the files above | derived | CC BY 4.0 (this project), with attribution to the underlying sources | – |
| `ukseat/config_2024.py` (poll numbers) | final 2024 polls: GB average of BPC members' final polls; Scotland and Wales means of each pollster's last poll | British Polling Council, "The performance of the polls in the 2024 general election"; pollsters' published tables as listed on Wikipedia's sub-national polling page | factual figures quoted with source | the BPC page exists and is dated; the sub-national inputs were not re-checked against each pollster's tables |
| `out/*`, `docs/*` | stored model outputs and the note | this project | CC BY 4.0 (attribution: Nari Han), subject to the source terms above | – |

**Attribution.**
- UK Parliament data (`data/gb2024_seats_v0_1.csv`, `data/old_bound_hist_v0_1.csv`, the nomination columns of
  `data/gb2024_v02_covariates.csv`, and the files derived from them): Contains Parliamentary information licensed under
  the Open Parliament Licence v3.0 (https://www.parliament.uk/site-information/copyright-parliament/open-parliament-licence/).
- BES 2024 constituency results file (the census and Leave columns of `data/gb2024_v02_covariates.csv`): Perrett, S.,
  Scott, R., Griffiths, J., Prosser, C., Bailey, J., Evans, G., Fieldhouse, E., Green, J. and Mellon, J. (2025). 2024
  BES Constituency Results with Census and Candidate Data. University of Manchester. doi:10.48420/28430672. Licensed
  under CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/). Changes: columns selected and renamed, age bands
  summed, merged with the seat list by ONS code (`scripts/prep_v02.py`).
- This project's note and stored outputs: CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/), attribution: Nari
  Han.

**Corrections.**
- v0.1 of the note described the BES file as "CC BY 4.0, doi:10.48420/26340568". That DOI belongs to another dataset: Miori (2024), notional 2014 Scottish independence referendum results, which is cited inside the BES codebook. The licence had not been checked at that time.
- v0.2 cited the version DOI `10.48420/28430672.v1`, which could not be confirmed. v0.2.1 cites the base DOI shown in the University of Manchester record.
