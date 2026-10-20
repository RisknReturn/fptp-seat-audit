"""v0.2 inputs known before polling day: 2021 census + Leave estimates (BES 2024 results file), 2024 nominations
(which parties stand; independents; sitting MPs standing as independents), strong 2019 'other' candidates re-standing."""
import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import sqlite3, pandas as pd, numpy as np
DB = _os.environ.get("UKSEAT_DB", _os.path.join(_ROOT, "..", "data_raw", "psephology.db"))
c = sqlite3.connect(DB)
seats = pd.read_csv(f"{_ROOT}/data/gb2024_seats_v0_1.csv", index_col=0)
bes = pd.read_stata(f"{_ROOT}/data/bes/BES-2024-General-Election-results-file-v1.0.dta")
cols = {"HanrettyLeave": "leave", "c21Muslim": "muslim", "c21Level4plus": "degree", "c21EthnicityWhiteBritish": "white_british"}
b = bes[["ONSConstID"] + list(cols)].rename(columns=cols).rename(columns={"ONSConstID": "ons"})
age = bes[["c21Age20to24", "c21Age25to29"]].sum(axis=1); b["age20_29"] = age.values
old = bes[["c21Age65to69", "c21Age70to74", "c21Age75to79", "c21Age80to84", "c21Age85plus"]].sum(axis=1); b["age65p"] = old.values
q = """select cg.constituency_area_id area_id, ca.candidate_given_name gn, ca.candidate_family_name fn, ca.candidate_is_sitting_mp mp,
       ca.is_standing_as_independent ind, ca.is_standing_as_commons_speaker spk, ce.political_party_id pid
       from candidacies ca join elections e on ca.election_id=e.id join constituency_groups cg on cg.id=e.constituency_group_id
       left join certifications ce on ce.candidacy_id=ca.id and ce.adjunct_to_certification_id is null where e.general_election_id={g}"""
c24 = pd.read_sql(q.format(g=6), c)
PM = {1: "Con", 3: "Lab", 4: "LD", 257: "Ref", 2: "Grn", 255: "Grn", 11: "SNP", 5: "PC"}
c24["party"] = c24.pid.map(PM).fillna("Oth")
stand = c24.pivot_table(index="area_id", columns="party", values="gn", aggfunc="count", fill_value=0)
out = seats[["name", "ons", "country", "region"]].copy()
for k in ["Con", "Lab", "LD", "Ref", "Grn", "SNP", "PC", "Oth"]:
    out["stand24_" + k] = (stand.reindex(out.index)[k].fillna(0) > 0) if k in stand else False
# sitting MPs standing as independents (nominations are public before polling day)
ind_mp = c24[(c24.mp == 1) & (c24.ind == 1)]
out["ind_sitting_mp"] = out.index.isin(ind_mp.area_id)
out["wpb24"] = out.index.isin(c24[c24.pid == 261].area_id)
# strong 2019 'other' (>=10%, GB, not Speaker) who stand again in 2024 (name match), mapped via overlaps
c19 = pd.read_sql(q.format(g=4) .replace("ca.is_standing_as_commons_speaker spk", "ca.is_standing_as_commons_speaker spk, ca.vote_share sh"), c)
c19["party"] = c19.pid.map({**PM, 153: "Ref", 8: "UKIP"}).fillna("Oth")
strong = c19[(c19.party == "Oth") & (c19.sh >= 0.10) & (c19.spk == 0)]
names24 = c24[c24.party == "Oth"][["area_id", "gn", "fn"]]
m = strong.merge(names24, on=["gn", "fn"], suffixes=("_19", "_24"))
out["oth_returning_share"] = 0.0
for r in m.itertuples():
    if r.area_id_24 in out.index:
        out.loc[r.area_id_24, "oth_returning_share"] = r.sh
print("returning strong others:", m[["gn", "fn", "area_id_24", "sh"]].values.tolist())
print("sitting MPs as independents:", ind_mp[["gn", "fn", "area_id"]].values.tolist())
out = out.reset_index().merge(b, on="ons", how="left").set_index("area_id")
print("missing census:", out.muslim.isna().sum(), "| stand counts:", {k: int(out['stand24_' + k].sum()) for k in ['Con','Lab','LD','Ref','Grn','SNP','PC','Oth']})
for v in ["leave", "muslim", "degree", "age20_29", "age65p", "white_british"]:
    out[v] = out[v].astype(float)
    print(v, out[v].describe()[["mean", "min", "max"]].round(3).to_dict())
out.to_csv(f"{_ROOT}/data/gb2024_v02_covariates.csv")
