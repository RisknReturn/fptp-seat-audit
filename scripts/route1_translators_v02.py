"""Route 1 (R2 follow-up): where does the 2024 Labour seat-efficiency gap come from?

All translators are given the ACTUAL 2024 vote totals (oracle; development diagnostic, not a forecast) and map the
2019 notional result of every seat (2024 boundaries) to 2024 seat shares. The largest 'other' candidate of each seat is
held at its actual 2024 share in every translator (independents are idiosyncratic; this isolates the party-system
translation). Seats are then counted (a) deterministically (largest share) and (b) with the v0.2 seat-level noise
(two-scale Gumbel mixture on log shares, the same win_probs as the engine, no regional factors).

Translators
  UNS_GB / UNS_country / UNS_region : additive (uniform) swing at GB / country / region level
  PS_country / PS_region            : proportional swing; Reform baseline = 2019 Brexit share, or the 2015 UKIP map
                                       where the Brexit Party did not stand (as in the v0.2 model)
  model_mu0 / model_mustar          : the v0.2.1 engine (oracle, regional factors integrated), from mu_inversion_v02
  ECO_country                       : ecological transfer matrix (2019 party -> 2024 party) per country, fitted to the
                                       2024 seat results (post-hoc reconstruction)
  ECO_country_x_context             : the same, separately for seat contexts (2019 winner x runner-up among Con/Lab/LD)
No turnout / abstention data are available in the seat file, so transfers are among voters only.
Writes out/route1_translators_v02.json.
"""
import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import sys, json, numpy as np, pandas as pd
from scipy.optimize import minimize
sys.path.insert(0, _ROOT)
src = open(f"{_ROOT}/scripts/v02_experiments.py").read()
sys.argv = [sys.argv[0], "4"]
exec(src[:src.index('out = {"kappa"')])

s = pd.read_csv(f"{_ROOT}/data/gb2024_seats_v0_1.csv", index_col=0)
s = s[(s.n19_Spk == 0) & (s.a24_Spk == 0)].copy()                   # drop the Speaker's seat
P = ["Con", "Lab", "LD", "Ref", "Grn", "SNP", "PC"]                  # parties translated


def shares(prefix):
    tot = s[[c for c in s.columns if c.startswith(prefix + "_") and not c.endswith("OthMax")]].sum(axis=1)
    out = pd.DataFrame({k: s[f"{prefix}_{k}"] for k in P}, index=s.index)
    out["Ref"] = out["Ref"] + s[f"{prefix}_UKIP"]
    out["OthMax"] = s[f"{prefix}_OthMax"]
    out["OthRest"] = s[f"{prefix}_Oth"] - s[f"{prefix}_OthMax"]
    return out.div(tot, axis=0), tot


v19, t19 = shares("n19")
v24, t24 = shares("a24")
oth24 = v24.OthMax.values                    # held at actual
party_mass24 = 1.0 - v24.OthMax.values - v24.OthRest.values
ACT = {k: int(v) for k, v in C.ACTUAL_SEATS_GB_2024.items()}
act_win = np.array(P + ["Oth"])[np.c_[v24[P].values, v24.OthMax.values].argmax(1)]
mtmp = mk(POLLS_ACT)
pi, tau1, tau2 = mtmp.pi, mtmp.tau1, mtmp.tau2


def seat_counts(sh):
    """sh: (n, 7) predicted party shares (each row summing to party_mass24). Returns deterministic and noisy counts."""
    full = np.c_[sh, oth24]
    win = np.array(P + ["Oth"])[full.argmax(1)]
    det = {k: int((win == k).sum()) for k in P + ["Oth"]}
    with np.errstate(divide="ignore"):
        eta = np.log(np.maximum(full, 0))
    pres = full > 0
    pr = 0.0
    for wgt, tau in ((1 - pi, tau1), (pi, tau2)):
        x = np.where(pres, eta / tau, -np.inf)
        mx = x.max(1, keepdims=True)
        e = np.exp(x - mx)
        pr = pr + wgt * e / e.sum(1, keepdims=True)
    exp_ = {k: float(pr[:, j].sum()) for j, k in enumerate(P + ["Oth"])}
    hit = float((win == act_win).mean())
    by_country = {c: {k: int(((win == k) & (s.country.values == c)).sum()) for k in ("Lab", "Con", "LD", "SNP", "Ref")}
                  for c in ("England", "Scotland", "Wales")}
    return dict(deterministic=det, deterministic_by_country=by_country, expected_with_seat_noise=exp_, hit_rate=hit,
                vote_mse=float(((sh - v24[P].values) ** 2).sum(1).mean()))


def renorm(x):
    x = np.maximum(x, 0)
    return x / x.sum(1, keepdims=True) * party_mass24[:, None]


def level_totals(v, t, groups):
    return (v[P].mul(t, axis=0)).groupby(groups).sum().div(t.groupby(groups).sum(), axis=0)


res = dict(note=__doc__, actual=ACT, n_seats=int(len(s)), translators={})
geo = {"GB": pd.Series("GB", index=s.index), "country": s.country, "region": s.region}
# --- uniform (additive) swing
for lvl in ("GB", "country", "region"):
    g = geo[lvl]
    d = level_totals(v24, t24, g) - level_totals(v19, t19, g)
    sh = v19[P].values + d.loc[g].values
    res["translators"][f"UNS_{lvl}"] = seat_counts(renorm(sh))
# --- proportional swing with the 2015 UKIP map as Reform baseline where the Brexit Party did not stand
base = v19[P].copy()
noref = base.Ref <= 0
u15 = s["s15map_UKIP"].fillna(0)
for ctry in s.country.unique():                    # Brexit-2019 / UKIP-2015 ratio where both are present, per country
    ok = (s.country == ctry) & (~noref) & (u15 > 0)
    kfac = float(np.median(base.Ref[ok] / u15[ok])) if ok.any() else 0.0
    sel = noref & (s.country == ctry)
    base.loc[sel, "Ref"] = kfac * u15[sel]
mass19 = v19[P].sum(1)
base = base.div(base.sum(1), axis=0).mul(mass19, axis=0)
for lvl in ("country", "region"):
    g = geo[lvl]
    B = level_totals(base.assign(), t19, g)
    ratio = (level_totals(v24, t24, g) / B).replace([np.inf, -np.inf], 0).fillna(0)
    sh = base[P].values * ratio.loc[g].values
    res["translators"][f"PS_{lvl}"] = seat_counts(renorm(sh))
# --- v0.2.1 engine (oracle), from the mu inversion run (regional factors integrated, rule (3,3,3,3,3,3))
mi = json.load(open(f"{_ROOT}/out/mu_inversion_v02.json"))
for tag, key in (("model_mu0", "mu=0"), ("model_h_obs", "mu=h_obs (v0.2 practice)"), ("model_mustar", "mu=mu_star (inverted)")):
    res["translators"][tag] = dict(expected_with_seat_noise=mi["oracle"][key]["mean_seats"],
                                   note="engine expected seats incl. regional factors; deterministic counts not computed")
# --- ecological transfer matrices (post-hoc)
SRC = P + ["OthAll"]
X19 = np.c_[v19[P].values, (v19.OthMax + v19.OthRest).values]       # 8 source shares (sum 1)
Y24 = v24[P].values / party_mass24[:, None]                           # 7 target shares among parties


def fit_T(ix, w):
    Xi, Yi, wi = X19[ix], Y24[ix], w[ix]
    present = Xi.sum(0) > 1e-9
    k_src, k_tgt = Xi.shape[1], Yi.shape[1]
    tgt_ok = Yi.sum(0) > 1e-9

    def unpack(th):
        Z = th.reshape(k_src, k_tgt)
        Z = np.where(tgt_ok[None, :], Z, -50.0)
        Z = Z - Z.max(1, keepdims=True)
        E = np.exp(Z)
        return E / E.sum(1, keepdims=True)

    def loss(th):
        T = unpack(th)
        pred = Xi @ T
        pred = pred / pred.sum(1, keepdims=True)
        return float((wi[:, None] * (pred - Yi) ** 2).sum() / wi.sum())

    th0 = np.zeros(k_src * k_tgt)
    for j in range(k_tgt):                                   # start near 'stay with the same party'
        if j < k_src:
            th0[j * k_tgt + j] = 2.0
    r = minimize(loss, th0, method="L-BFGS-B", options=dict(maxiter=5000))
    T = unpack(r.x)
    T[~present] = np.nan
    return T, r


def eco(groups, tag):
    pred = np.zeros_like(Y24)
    mats = {}
    w = t24.values.astype(float)
    for gname in sorted(set(groups)):
        ix = np.where(groups == gname)[0]
        T, r = fit_T(ix, w)
        Tz = np.nan_to_num(T)
        p = X19[ix] @ Tz
        pred[ix] = p / p.sum(1, keepdims=True)
        mats[gname] = dict(n=int(len(ix)), loss=float(r.fun), success=bool(r.success),
                           T=pd.DataFrame(T, index=SRC, columns=P).round(4).to_dict())
    out = seat_counts(pred * party_mass24[:, None])
    out["matrices"] = mats
    res["translators"][tag] = out


eco(s.country.values, "ECO_country")
# seat context = 2019 winner > runner-up. England/Wales among Con/Lab/LD; Scotland among SNP/Lab/Con/LD.
# Groups with fewer than 20 seats are merged into '<country>:other' (Wales is one group: 32 seats).
def context_labels():
    lab = []
    for i, c in enumerate(s.country.values):
        if c == "Wales":
            lab.append("Wales"); continue
        cand = ["SNP", "Lab", "Con", "LD"] if c == "Scotland" else ["Con", "Lab", "LD"]
        vals = v19.iloc[i][cand].values
        o = np.argsort(-vals)
        lab.append(f"{c}:{cand[o[0]]}>{cand[o[1]]}")
    lab = pd.Series(lab, index=s.index)
    small = lab.map(lab.value_counts()) < 20
    lab[small] = s.country[small] + ":other"
    return lab.values


ctx = context_labels()
res["context_sizes"] = pd.Series(ctx).value_counts().to_dict()
eco(ctx, "ECO_country_x_context")

res["actual_by_country"] = {c: {k: int(((act_win == k) & (s.country.values == c)).sum()) for k in ("Lab", "Con", "LD", "SNP", "Ref")}
                            for c in ("England", "Scotland", "Wales")}
print(f"{'translator':24s} {'Lab det':>8s} {'Lab exp':>8s} {'Con':>6s} {'LD':>6s} {'Ref':>6s} {'SNP':>6s} {'hit':>6s} {'voteMSE':>8s}")
for tag, r in res["translators"].items():
    d = r.get("deterministic", {}); e = r["expected_with_seat_noise"]
    print(f"{tag:24s} {d.get('Lab', float('nan')):>8} {e['Lab']:8.1f} {e['Con']:6.1f} {e['LD']:6.1f} {e['Ref']:6.1f} {e['SNP']:6.1f} "
          f"{r.get('hit_rate', float('nan')):6.3f} {r.get('vote_mse', float('nan')):8.5f}", flush=True)
print("actual", ACT)
for tag in ("UNS_country", "PS_country", "ECO_country", "ECO_country_x_context"):
    print(tag, res["translators"][tag]["deterministic_by_country"], flush=True)
print("actual by country", res["actual_by_country"])
json.dump(res, open(f"{_ROOT}/out/route1_translators_v02.json", "w"), indent=1, default=float)

# --- out-of-sample check: K-fold by seat, stratified by group (fit T on the other folds, predict held-out seats)
def eco_cv(groups, tag, K=5, seed=20261020):
    rng = np.random.default_rng(seed)
    w = t24.values.astype(float)
    fold = np.empty(len(groups), int)
    for gname in sorted(set(groups)):          # sorted: fold assignment must not depend on string hashing
        ix = np.where(groups == gname)[0]
        fold[ix] = rng.permutation(np.arange(len(ix)) % K)
    pred = np.zeros_like(Y24)
    for f in range(K):
        for gname in sorted(set(groups)):
            tr = np.where((groups == gname) & (fold != f))[0]
            te = np.where((groups == gname) & (fold == f))[0]
            if len(te) == 0:
                continue
            T, _ = fit_T(tr, w)
            p = X19[te] @ np.nan_to_num(T)
            pred[te] = p / np.maximum(p.sum(1, keepdims=True), 1e-12)
    out = seat_counts(pred * party_mass24[:, None])
    out["note"] = f"{K}-fold out-of-sample by seat, stratified by group, seed {seed}"
    res["translators"][tag] = out


eco_cv(s.country.values, "ECO_country_CV5")
eco_cv(ctx, "ECO_country_x_context_CV5")
for tag in ("ECO_country_CV5", "ECO_country_x_context_CV5"):
    r = res["translators"][tag]; e = r["expected_with_seat_noise"]
    print(f"{tag:26s} Lab det {r['deterministic']['Lab']} exp {e['Lab']:.1f} Con {e['Con']:.1f} LD {e['LD']:.1f} SNP {e['SNP']:.1f} "
          f"hit {r['hit_rate']:.3f} voteMSE {r['vote_mse']:.5f} | {r['deterministic_by_country']}", flush=True)
print("contexts", res["context_sizes"])
json.dump(res, open(f"{_ROOT}/out/route1_translators_v02.json", "w"), indent=1, default=float)

# --- robustness: (i) contexts by 2019 winner only, (ii) context matrices shrunk toward the country matrix
def fit_T_prior(ix, w, T0, alpha):
    """As fit_T, plus alpha * ||T - T0||^2 (T0 = the country matrix; rows absent in T0 are not penalised)."""
    Xi, Yi, wi = X19[ix], Y24[ix], w[ix]
    present = Xi.sum(0) > 1e-9
    k_src, k_tgt = Xi.shape[1], Yi.shape[1]
    tgt_ok = Yi.sum(0) > 1e-9
    T0z = np.nan_to_num(T0)
    mask = ~np.isnan(T0).any(1)

    def unpack(th):
        Z = th.reshape(k_src, k_tgt)
        Z = np.where(tgt_ok[None, :], Z, -50.0)
        Z = Z - Z.max(1, keepdims=True)
        E = np.exp(Z)
        return E / E.sum(1, keepdims=True)

    def loss(th):
        T = unpack(th)
        pred = Xi @ T
        pred = pred / pred.sum(1, keepdims=True)
        fit = (wi[:, None] * (pred - Yi) ** 2).sum() / wi.sum()
        return float(fit + alpha * ((T - T0z)[mask & present] ** 2).sum())

    th0 = np.log(np.maximum(T0z, 1e-6)).ravel()
    r = minimize(loss, th0, method="L-BFGS-B", options=dict(maxiter=5000))
    T = unpack(r.x); T[~present] = np.nan
    return T


def eco_cv_prior(groups, tag, alpha, K=5, seed=20261020):
    rng = np.random.default_rng(seed)
    w = t24.values.astype(float)
    fold = np.empty(len(groups), int)
    for gname in sorted(set(groups)):          # sorted: fold assignment must not depend on string hashing
        ix = np.where(groups == gname)[0]
        fold[ix] = rng.permutation(np.arange(len(ix)) % K)
    pred = np.zeros_like(Y24)
    for f in range(K):
        for ctry in ("England", "Scotland", "Wales"):
            trc = np.where((s.country.values == ctry) & (fold != f))[0]
            T0, _ = fit_T(trc, w)
            for gname in sorted(set(groups[s.country.values == ctry])):
                tr = np.where((groups == gname) & (fold != f))[0]
                te = np.where((groups == gname) & (fold == f))[0]
                if len(te) == 0:
                    continue
                T = fit_T_prior(tr, w, T0, alpha) if alpha < np.inf else T0
                p = X19[te] @ np.nan_to_num(T)
                pred[te] = p / np.maximum(p.sum(1, keepdims=True), 1e-12)
    out = seat_counts(pred * party_mass24[:, None])
    out["note"] = f"{K}-fold out-of-sample, context matrices shrunk toward the country matrix, alpha={alpha}"
    res["translators"][tag] = out


win_ctx = np.array([f"{c}:{w}" for c, w in zip(s.country.values,
                    np.where(s.country.values == "Scotland",
                             np.array(["SNP", "Lab", "Con", "LD"])[v19[["SNP", "Lab", "Con", "LD"]].values.argmax(1)],
                             np.array(["Con", "Lab", "LD"])[v19[["Con", "Lab", "LD"]].values.argmax(1)]))])
wl = pd.Series(win_ctx)
win_ctx = np.where(wl.map(wl.value_counts()) < 20, s.country.values + ":other", win_ctx)
eco_cv(win_ctx, "ECO_country_x_winner_CV5")
for a_ in (0.003, 0.03, 0.3):
    eco_cv_prior(ctx, f"ECO_context_shrunk_a{a_}_CV5", a_)
for tag in ["ECO_country_x_winner_CV5"] + [f"ECO_context_shrunk_a{a_}_CV5" for a_ in (0.003, 0.03, 0.3)]:
    r = res["translators"][tag]; e = r["expected_with_seat_noise"]
    print(f"{tag:30s} Lab det {r['deterministic']['Lab']} exp {e['Lab']:.1f} Con {e['Con']:.1f} LD {e['LD']:.1f} SNP {e['SNP']:.1f} "
          f"hit {r['hit_rate']:.3f} | Eng Lab {r['deterministic_by_country']['England']['Lab']} Sco Lab {r['deterministic_by_country']['Scotland']['Lab']}", flush=True)
json.dump(res, open(f"{_ROOT}/out/route1_translators_v02.json", "w"), indent=1, default=float)

# --- combination: England/Wales by winner x runner-up, Scotland by 2019 winner
ctx_mix = np.where(s.country.values == "Scotland", win_ctx, ctx)
eco_cv(ctx_mix, "ECO_mixed_context_CV5")
r = res["translators"]["ECO_mixed_context_CV5"]; e = r["expected_with_seat_noise"]
print(f"ECO_mixed_context_CV5 Lab det {r['deterministic']['Lab']} exp {e['Lab']:.1f} Con {e['Con']:.1f} LD {e['LD']:.1f} SNP {e['SNP']:.1f} "
      f"hit {r['hit_rate']:.3f} | {r['deterministic_by_country']}", flush=True)
json.dump(res, open(f"{_ROOT}/out/route1_translators_v02.json", "w"), indent=1, default=float)


# --- fold-assignment sensitivity of the main variants (three seeds)
res["fold_seed_sensitivity"] = {}
for seed_ in (1, 2, 3):
    for g_, tag_ in ((ctx_mix, "mixed"), (ctx, "context"), (s.country.values, "country")):
        key = f"{tag_}_seed{seed_}"
        eco_cv(g_, "tmp", seed=seed_)
        t_ = res["translators"].pop("tmp")
        res["fold_seed_sensitivity"][key] = dict(Lab_det=t_["deterministic"]["Lab"], Lab_exp=t_["expected_with_seat_noise"]["Lab"],
                                                 Con_exp=t_["expected_with_seat_noise"]["Con"], LD_exp=t_["expected_with_seat_noise"]["LD"],
                                                 SNP_exp=t_["expected_with_seat_noise"]["SNP"], hit=t_["hit_rate"])
print("fold seeds", json.dumps(res["fold_seed_sensitivity"]), flush=True)
json.dump(res, open(f"{_ROOT}/out/route1_translators_v02.json", "w"), indent=1, default=float)
