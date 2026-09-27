"""Regenerate the offline catalog from a pinned, licensed Biopython release.

Parse literals only: downloaded Python is never executed. No sequence data is sent.
"""
import ast
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen

BASE = "https://raw.githubusercontent.com/biopython/biopython/biopython-187/"
ROOT = Path(__file__).resolve().parents[1]


def main():
    url = BASE + "Bio/Restriction/Restriction_Dictionary.py"
    raw = urlopen(url, timeout=30).read()
    records = []
    for node in ast.parse(raw.decode("utf-8")).body:
        if not isinstance(node, ast.Assign) or len(node.targets) != 1:
            continue
        target = node.targets[0]
        if not isinstance(target, ast.Subscript) or not isinstance(target.value, ast.Name) or target.value.id != "rest_dict":
            continue
        name = ast.literal_eval(target.slice)
        entry = ast.literal_eval(node.value)
        size = len(entry["site"])
        cuts = [[entry["fst5"], size + entry["fst3"]]] if entry["fst5"] is not None and entry["fst3"] is not None else []
        if entry["scd5"] is not None and entry["scd3"] is not None:
            cuts.append([entry["scd5"], size + entry["scd3"]])
        records.append({"name": name, "recognition": entry["site"], "cuts": cuts,
                        "suppliers": list(entry["suppl"]), "rebase_id": entry["id"]})
    neb_url = "https://www.neb.com/en-ca/tools-and-resources/selection-charts/isoschizomers"
    # Explicitly reviewed product names; do not infer arbitrary suffixes.
    products = {
        "PsiI": ["PsiI-v2"], "KpnI": ["KpnI-HF"], "ApoI": ["ApoI-HF"],
        "DraIII": ["DraIII-HF"], "AgeI": ["AgeI-HF"], "SpeI": ["SpeI-HF"],
        "AleI": ["AleI-v2"],
        "EcoRI": ["EcoRI-HF"], "BamHI": ["BamHI-HF"],
        "HindIII": ["HindIII-HF"], "BsaI": ["BsaI-HF", "BsaI-HFv2"],
    }
    known = {entry["name"]: entry for entry in records}
    for base, names in products.items():
        if base in known:
            known[base]["product_names"] = names
    data = {"source": url, "source_sha256": hashlib.sha256(raw).hexdigest(),
            "product_names_source": neb_url, "product_names_checked": "2026-09-27",
            "additional_product_names_source": "https://enzymefinder.neb.com/",
            "release": "Biopython 1.87", "rebase_version": "404 (2024)",
            "enzymes": records}
    out = ROOT / "primerblast_oss" / "restriction_catalog.json"
    out.write_text(json.dumps(data, ensure_ascii=False, separators=(",", ":")) + "\n", encoding="utf-8")
    notice = ROOT / "THIRD_PARTY_BIOPYTHON_LICENSE.rst"
    notice.write_bytes(b"Restriction dictionary data derived from Biopython.\nCopyright (C) 2004. Frederic Sohm.\nSource: " + url.encode("utf-8") + b"\n\n" + urlopen(BASE + "LICENSE.rst", timeout=30).read())
    print("Catalog records:", len(records))


if __name__ == "__main__":
    main()
