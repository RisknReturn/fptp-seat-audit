"""Inputs for the 2024 backtest, all available before polling day (4 July 2024),
except the 'actual' block used only for evaluation."""
import numpy as np

# GB final polls: British Polling Council, 'The performance of the polls in the 2024 general election'
# (average of 19 BPC members' final polls)
POLLS_GB_2024 = {"Con": 21.4, "Lab": 38.6, "LD": 11.3, "Grn": 6.5, "Ref": 16.3, "Oth": 6.0}
# Scotland: mean of the last poll of each pollster, fieldwork ending 20 Jun - 2 Jul 2024
# (Savanta, More in Common, Redfield & Wilton, Norstat, Survation, YouGov; Wikipedia sub-national polling page)
_S = np.array([[34, 31, 15, 9, 3, 6], [30, 35, 16, 9, 2, 7], [32, 38, 11, 7, 3, 8],
               [31, 35, 13, 8, 3, 8], [31, 37, 14, 7, 3, 8], [29, 35, 11, 11, 5, 8]], float).mean(0)
POLLS_SCOT_2024 = dict(zip(["SNP", "Lab", "Con", "LD", "Grn", "Ref"], _S))
# Wales: last poll of each pollster, fieldwork ending 14 Jun - 1 Jul 2024 (YouGov, More in Common, R&W, Savanta)
_W = np.array([[40, 16, 14, 7, 5, 16], [42, 22, 9, 4, 5, 14], [46, 15, 10, 7, 4, 17], [49, 19, 12, 5, 3, 12]], float).mean(0)
POLLS_WALES_2024 = dict(zip(["Lab", "Con", "PC", "LD", "Grn", "Ref"], _W))
POLLS_2024 = {"GB": POLLS_GB_2024, "Scotland": POLLS_SCOT_2024, "Wales": POLLS_WALES_2024}

# Final-poll-average errors (poll minus GB result, points), elections before 2024.
# 2010: Wikipedia final-poll table (mean of 9 pollsters' last polls). 2015: mean of the 10 final polls of
# BPC members (Wikipedia table; excludes SurveyMonkey and Survation's unpublished 6 May poll), consistent
# with the Sturgis inquiry's "34% / 34%". 2017, 2019: BPC final-poll tables.
POLL_ERR = {
    2010: {"Con": 35.56 - 36.9, "Lab": 27.56 - 29.7, "LD": 27.44 - 23.6},
    2015: {"Con": 33.7 - 37.8, "Lab": 33.7 - 31.2, "LD": 9.1 - 8.1, "Ref": 12.6 - 12.9, "Grn": 4.8 - 3.8},
    2017: {"Con": 43.7 - 43.5, "Lab": 35.8 - 41.0, "LD": 7.9 - 7.6, "Ref": 4.4 - 1.9, "Grn": 2.1 - 1.6},
    2019: {"Con": 43.3 - 44.7, "Lab": 33.5 - 33.0, "LD": 11.7 - 11.8, "Ref": 2.9 - 2.1, "Grn": 2.9 - 2.8},
}
# national polling-error sd (points): RMS of 2010-2019 errors for Con/Lab/LD; floors for Ref/Grn
# because their history is dominated by standing-pattern artefacts (see design doc).
SIGMA_NAT = {"Con": 2.3, "Lab": 3.0, "LD": 2.0, "Ref": 2.0, "Grn": 1.0}
# Correlation: start from corr(Con,Lab) = -0.3 (2010-19 errors) with 'Other' sd 1.0, then condition on
# the errors summing to zero across all six groups (shares add to 100), and keep the conditional
# correlation with the marginal sds above. Implied 'Other' error sd = 1.1 points.
def _corr_sum_zero(sd5, sd_oth=1.0, rho_cl=-0.3):
    sd = np.array(list(sd5) + [sd_oth]); R = np.eye(6); R[0, 1] = R[1, 0] = rho_cl
    S = np.outer(sd, sd) * R; one = np.ones(6)
    Sc = S - np.outer(S @ one, one @ S) / (one @ S @ one)
    d = np.sqrt(np.diag(Sc))
    return (Sc / np.outer(d, d))[:5, :5]
CORR_NAT = _corr_sum_zero([SIGMA_NAT[k] for k in ["Con", "Lab", "LD", "Ref", "Grn"]])
# bias reference vector for lambda scenarios: the 2015 miss ("shy Tory"), poll minus result
E_REF_2015 = POLL_ERR[2015]

ACTUAL_GB_2024 = {"Con": 24.4, "Lab": 34.7, "LD": 12.5, "Grn": 6.9, "Ref": 14.7, "Oth": 6.8}
ACTUAL_SEATS_GB_2024 = {"Lab": 411, "Con": 121, "LD": 72, "SNP": 9, "Ref": 5, "Oth": 5, "Grn": 4, "PC": 4, "Spk": 1}
