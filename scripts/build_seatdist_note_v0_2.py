"""Assemble docs/SEATDIST_NOTE_v0_2.md from the template and the generated tables (out/seatdist_note_tables_v0_2.md).
Placeholders {{M1}}, {{N1}} ... {{N9}} in the template are replaced by the matching '### Table ...' blocks."""
import os, re
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
tab = open(f"{ROOT}/out/seatdist_note_tables_v0_2.md").read()
blocks = {}
for m in re.finditer(r"^### Table ([MN]\d+[abc]?)\.(.*?)(?=^### Table |\Z)", tab, re.S | re.M):
    blocks[m.group(1)] = ("#### Table " + m.group(1) + "." + m.group(2)).rstrip() + "\n"
tpl = open(f"{ROOT}/docs/SEATDIST_NOTE_v0_2.template.md").read()
need = set(re.findall(r"\{\{([MN]\d+[abc]?)\}\}", tpl))
missing = need - set(blocks)
assert not missing, missing
out = re.sub(r"\{\{([MN]\d+[abc]?)\}\}", lambda m: blocks[m.group(1)], tpl)
open(f"{ROOT}/docs/SEATDIST_NOTE_v0_2.md", "w").write(out)
print("tables used:", sorted(need), "| unused:", sorted(set(blocks) - need), "| words:", len(out.split()))
