"""Assemble docs/SEATDIST_NOTE_v0_2_2.md from the template and the generated tables (out/seatdist_note_tables_v0_2.md
for M1 and N1-N9, out/seatdist_note_tables_v0_2_2.md for N7d). Placeholders {{M1}}, {{N1}} ... in the template are
replaced by the matching '### Table ...' blocks."""
import os, re
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
blocks = {}
for f in ("seatdist_note_tables_v0_2.md", "seatdist_note_tables_v0_2_2.md"):
    tab = open(f"{ROOT}/out/{f}").read()
    for m in re.finditer(r"^### Table ([MN]\d+[a-d]?)\.(.*?)(?=^### Table |\Z)", tab, re.S | re.M):
        assert m.group(1) not in blocks, m.group(1)
        blocks[m.group(1)] = ("#### Table " + m.group(1) + "." + m.group(2)).rstrip() + "\n"
tpl = open(f"{ROOT}/docs/SEATDIST_NOTE_v0_2_2.template.md").read()
need = set(re.findall(r"\{\{([MN]\d+[a-d]?)\}\}", tpl))
missing = need - set(blocks)
assert not missing, missing
out = re.sub(r"\{\{([MN]\d+[a-d]?)\}\}", lambda m: blocks[m.group(1)], tpl)
assert "{{" not in out, re.findall(r"\{\{[^}]*\}\}", out)
open(f"{ROOT}/docs/SEATDIST_NOTE_v0_2_2.md", "w").write(out)
print("tables used:", sorted(need), "| unused:", sorted(set(blocks) - need), "| words:", len(out.split()))
