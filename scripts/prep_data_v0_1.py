"""Build constituency-level tables from the UK Parliament psephology database.

Source: https://github.com/ukparliament/psephology-datasette (psephology.db), the
House of Commons Library election results data behind electionresults.parliament.uk.
Outputs (data/):
  gb2024_seats_v0_1.csv   : GB seats on 2024 boundaries, 2019 notional + 2024 actual votes,
                             UKIP 2015 share mapped to 2024 boundaries.
  old_bound_hist_v0_1.csv : 2010/2015/2017/2019 actual votes on 2010 boundaries (calibration).
"""
import os as _os
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import sqlite3, sys, pandas as pd, numpy as np

DB = sys.argv[1] if len(sys.argv) > 1 else _os.environ.get("UKSEAT_DB", _os.path.join(_ROOT, "..", "data_raw", "psephology.db"))
OUT = f"{_ROOT}/data"
c = sqlite3.connect(DB)

PARTY_MAP = {1: "Con", 3: "Lab", 4: "LD", 257: "Ref", 153: "Ref", 8: "UKIP",
             2: "Grn", 255: "Grn", 267: "Grn", 11: "SNP", 5: "PC"}
COLS = ["Con", "Lab", "LD", "Ref", "UKIP", "Grn", "SNP", "PC", "Spk", "Oth"]
WCOLS = ["Con", "Lab", "LD", "Ref", "Grn", "SNP", "PC", "Spk", "OthMax"]

def votes(ge_id):
    q = f"""select e.id election_id, cg.constituency_area_id area_id, ca.vote_count v,
                   ca.is_standing_as_commons_speaker spk, ce.political_party_id pid
            from candidacies ca join elections e on ca.election_id=e.id
            join constituency_groups cg on cg.id=e.constituency_group_id
            left join certifications ce on ce.candidacy_id=ca.id and ce.adjunct_to_certification_id is null
            where e.general_election_id={ge_id}"""
    d = pd.read_sql(q, c)
    d["party"] = d["pid"].map(PARTY_MAP).fillna("Oth")
    d.loc[d.spk == 1, "party"] = "Spk"
    t = d.pivot_table(index="area_id", columns="party", values="v", aggfunc="sum", fill_value=0)
    for k in COLS:
        if k not in t:
            t[k] = 0
    t = t[COLS].astype(float)
    # largest single "other" candidate: an aggregate of several others cannot win a seat
    om = d[d.party == "Oth"].groupby("area_id").v.max()
    t["OthMax"] = om.reindex(t.index).fillna(0).astype(float)
    return t

areas = pd.read_sql("""select a.id area_id, a.name, a.geographic_code ons, a.boundary_set_id bset,
                              co.name country, r.name region
                       from constituency_areas a join countries co on co.id=a.country_id
                       left join english_regions r on r.id=a.english_region_id""", c).set_index("area_id")
areas["region"] = areas["region"].fillna(areas["country"])

# ---- 2024 boundaries: 2019 notional (GE 5) and 2024 actual (GE 6)
n19 = votes(5)
a24 = votes(6)
assert set(n19.index) == set(a24.index), "notional/actual area mismatch"
new = areas.loc[a24.index]
gb = new.country != "Northern Ireland"
print("2024-boundary seats:", len(new), "GB:", int(gb.sum()))

# ---- UKIP 2015 share mapped to 2024 boundaries through population overlaps
ov = pd.read_sql("select from_constituency_area_id f, to_constituency_area_id t, to_constituency_population tp from constituency_area_overlaps", c)
v15 = votes(2)
sh15 = v15[COLS].div(v15[COLS].sum(axis=1), axis=0)
ov = ov[ov.f.isin(sh15.index) & ov.t.isin(new.index)]
cov = ov.groupby("t").tp.sum()
print("overlap weight per new seat: min %.3f max %.3f, n=%d" % (cov.min(), cov.max(), len(cov)))
m = ov.merge(sh15, left_on="f", right_index=True)
mapped = {}
for k in ["UKIP", "Grn", "LD", "Con", "Lab"]:
    mapped[k] = (m.tp * m[k]).groupby(m.t).sum() / cov
mapped = pd.DataFrame(mapped)

out = new[["name", "ons", "country", "region"]].copy()
for k in COLS + ["OthMax"]:
    out["n19_" + k] = n19[k]
for k in COLS + ["OthMax"]:
    out["a24_" + k] = a24[k]
for k in mapped:
    out["s15map_" + k] = mapped[k]
out = out[gb].sort_values(["country", "region", "name"])
out.index.name = "area_id"
out.to_csv(f"{OUT}/gb2024_seats_v0_1.csv")
print("written", len(out), "GB seats; missing UKIP15 map:", int(out.s15map_UKIP.isna().sum()))

# ---- 2010-boundary history for calibration (GE 1..4)
hist = []
for ge, yr in [(1, 2010), (2, 2015), (3, 2017), (4, 2019)]:
    t = votes(ge)
    t.columns = [f"{yr}_{k}" for k in t.columns]
    hist.append(t)
h = pd.concat(hist, axis=1)
h = areas[["name", "country", "region"]].join(h, how="inner")
h = h[h.country != "Northern Ireland"]
h.to_csv(f"{OUT}/old_bound_hist_v0_1.csv")
print("history seats:", len(h))

# ---- sanity totals
for tag, t in [("2019 notional GB", out.filter(like="n19_")), ("2024 GB", out.filter(like="a24_"))]:
    s = t.sum(); s = s / s.sum() * 100
    print(tag, " ".join(f"{k.split('_')[1]} {v:.1f}" for k, v in s.items()))
for tag in ["n19", "a24"]:
    w = out[[f"{tag}_{k}" for k in WCOLS]].idxmax(axis=1).str[4:].value_counts()
    print(tag, "GB winners:", {k: int(v) for k, v in w.items()})
