"""Fill the Zenodo DOIs and the release date into the v0.2.3 files, regenerate docs/SEATDIST_NOTE_v0_2_3.md, and check
that no placeholder is left anywhere in the repository.

    python scripts/fill_doi_v0_2_3.py --concept 10.5281/zenodo.N --v022 10.5281/zenodo.M --date YYYY-MM-DD

--concept is the concept DOI (all versions), --v022 the version DOI of release v0.2.2 (v0.2.3 cannot carry its own DOI)
and --date the release date of v0.2.3. Placeholders: __CONCEPT_DOI__, __V022_DOI__, __V023_DATE__. CITATION.cff is left
with exactly one top-level date-released line (the v0.2.3 date), also after a merge that brought in the v0.2.2 date;
files with unresolved merge conflicts are refused before anything is written. Running it again with the same values
changes nothing; with other values it reports failure, because the files then carry DOIs that differ from the
arguments (restore them with `git checkout -- .`). Exit status 0 only if every check passes."""
import os, re, sys, json, shutil, argparse, datetime, subprocess
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FILES = ["CITATION.cff", "README.md", "README_ko.md", "docs/SEATDIST_NOTE_v0_2_3.template.md"]
NOTE = "docs/SEATDIST_NOTE_v0_2_3.md"
NEED = {"--concept": FILES + [NOTE], "--v022": FILES + [NOTE], "--date": ["CITATION.cff"]}  # where each value must be
DOI = re.compile(r"^10\.\d{4,9}/[^\s\"'<>()\[\]]+$")
SKIP_DIRS = {".git", "__pycache__", ".venv", "venv", ".ipynb_checkpoints"}

ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
ap.add_argument("--concept", required=True, help="concept DOI (all versions), e.g. 10.5281/zenodo.1234567")
ap.add_argument("--v022", required=True, help="version DOI of release v0.2.2")
ap.add_argument("--date", required=True, help="release date of v0.2.3, YYYY-MM-DD")
a = ap.parse_args()
for flag, v in (("--concept", a.concept), ("--v022", a.v022)):
    if not DOI.match(v):
        sys.exit(f"{flag}: expected a bare DOI such as 10.5281/zenodo.1234567, got {v!r}")
if a.concept == a.v022:
    sys.exit("--concept and --v022 must differ (concept DOI versus the v0.2.2 version DOI)")
try:
    assert re.fullmatch(r"\d{4}-\d{2}-\d{2}", a.date)
    datetime.date.fromisoformat(a.date)
except (AssertionError, ValueError):
    sys.exit(f"--date: expected YYYY-MM-DD, got {a.date!r}")
SUB = {"__CONCEPT_DOI__": a.concept, "__V022_DOI__": a.v022, "__V023_DATE__": a.date}
VAL = {"--concept": a.concept, "--v022": a.v022, "--date": a.date}


def read(rel):
    with open(f"{ROOT}/{rel}", encoding="utf-8", newline="") as f:  # newline="" keeps the files' own line endings
        return f.read()


def write(rel, text):
    with open(f"{ROOT}/{rel}", "w", encoding="utf-8", newline="") as f:
        f.write(text)


def one_date_released(cff):
    """Keep exactly one top-level date-released line, the v0.2.3 date, directly after the version line."""
    nl = "\r\n" if "\r\n" in cff else "\n"
    lines = [l for l in cff.split(nl) if not l.startswith("date-released:")]
    iv = [i for i, l in enumerate(lines) if l.startswith("version:")]
    assert len(iv) == 1, f"CITATION.cff: expected one top-level version line, found {len(iv)}"
    lines.insert(iv[0] + 1, f'date-released: "{a.date}"')
    return nl.join(lines)


# 0. refuse to touch files with unresolved merge conflicts (e.g. CITATION.cff after merging the v0.2.2 date line)
for rel in FILES:
    if re.search(r"^(<{7}|={7}|>{7})( |$)", read(rel), re.M):
        sys.exit(f"{rel}: unresolved merge conflict markers; resolve them first (keep this branch's lines), then rerun")

# 1. substitute
for rel in FILES:
    t = read(rel)
    n = {k: t.count(k) for k in SUB}
    for k, v in SUB.items():
        t = t.replace(k, v)
    if rel == "CITATION.cff":
        t = one_date_released(t)
    write(rel, t)
    print(f"{rel}: replaced {n}")

# 2. regenerate the v0.2.3 note from its template
subprocess.run([sys.executable, f"{ROOT}/scripts/build_seatdist_note_v0_2_3.py"], check=True)

# 3. checks
bad = []
for flag, rels in NEED.items():
    for rel in rels:
        if VAL[flag] not in read(rel):
            bad.append(f"{rel}: {flag} value {VAL[flag]} not found (filled earlier with other values?)")
cff = read("CITATION.cff")
top = lambda key: re.findall(rf"^{key}:\s*(.*?)\s*$", cff, re.M)
if top("version") != ["0.2.3"]:
    bad.append(f"CITATION.cff: top-level version is {top('version')}, expected ['0.2.3']")
if top("date-released") != [f'"{a.date}"']:
    bad.append(f"CITATION.cff: top-level date-released is {top('date-released')}, expected one line with {a.date}")
if top("doi") != [a.concept]:
    bad.append(f"CITATION.cff: top-level doi is {top('doi')}, expected {a.concept}")
zv = json.loads(read(".zenodo.json")).get("version")
if zv != "0.2.3":
    bad.append(f".zenodo.json: version is {zv!r}, expected '0.2.3'")
me = os.path.abspath(__file__)
for d, dirs, files in os.walk(ROOT):
    dirs[:] = [x for x in dirs if x not in SKIP_DIRS]
    for fn in files:
        p = os.path.join(d, fn)
        if os.path.abspath(p) == me:
            continue
        with open(p, "rb") as f:
            b = f.read()
        for k in SUB:
            for m in re.finditer(re.escape(k.encode()), b):
                line = b.count(b"\n", 0, m.start()) + 1
                bad.append(f"{os.path.relpath(p, ROOT)}:{line}: placeholder {k} left")
if shutil.which("cffconvert"):
    r = subprocess.run(["cffconvert", "--validate", "-i", f"{ROOT}/CITATION.cff"], capture_output=True, text=True)
    print("cffconvert:", (r.stdout + r.stderr).strip() or f"exit {r.returncode}")
    if r.returncode != 0:
        bad.append("CITATION.cff: cffconvert --validate failed")
else:
    print("cffconvert not installed; run `cffconvert --validate -i CITATION.cff` before committing")
if bad:
    print("FAILED", *bad, sep="\n  ")
    sys.exit(1)
print("OK: no placeholder left; review `git diff`, run the CI steps, then commit")
