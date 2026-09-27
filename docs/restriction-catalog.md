# Restriction catalog / 制限酵素カタログ

CAPS = Cleaved Amplified Polymorphic Sequence: genotype discrimination by
restriction digestion of PCR products (a PCR-RFLP approach). dCAPS = derived
CAPS: an intentional primer mismatch creates a diagnostic restriction site.
These are in-silico predictions until validated experimentally.

The offline catalog contains **1,088 enzyme names**, including noncommercial
entries, from Biopython 1.87's REBASE EMBOSS **404 (2024)** snapshot. This is a
broad catalog, not a claim to contain every enzyme discovered subsequently or
every vendor product. Supplier availability is historical, not current stock.
Verified NEB product variants are additional searchable names; arbitrary suffixes
are not silently removed.

- Same recognition **and cleavage**: related enzymes shown separately from product names.
- Same recognition, **different cleavage** (neoschizomers): never merged for prediction.
- Unknown cleavage: no fabricated cut diagram or fragment prediction.
- Two-pair cleavage: both pairs displayed in the catalog; excluded from the current single-pair CAPS/dCAPS model.
- Methylation-dependent/special substrates: excluded from ordinary PCR automatic candidates.

Default CAPS/dCAPS scans now consider **244 distinct cleavage geometries**, from
commercial-in-snapshot, supported enzymes and the existing curated panel.
Equivalent enzymes share computation and appear as related names in results.
This does not establish interchangeable methylation sensitivity, temperature,
buffer, multiple-site requirements, efficiency or suitability at DNA ends.
Always consult the selected product's documentation; wet validation remains open.

切断図の `|` は両鎖の切断境界です。認識配列外の N は隣接配列を表します。
切断座標は認識配列先頭からの0始まり境界で、配列外では負数や配列長を
超える値になります。上鎖と下鎖を同じ向きに並べ、5′・3′突出／平滑末端を
表示します。結果ではこれに加えて実際のアレル配列と産物内の切断位置を示します。

All catalog reads, searches and scans are local. No genome is uploaded to obtain
enzyme information. Runtime remains standard-library-only. The JSON is included
in installed packages; third-party license is retained.

Regenerate the public snapshot with `python tools/update_restriction_catalog.py`.
It downloads a pinned source and license and parses literal data without executing
downloaded Python. Provenance and SHA-256 are stored in the JSON. Product-name
additions require review against the cited primary sources.

Sources:

- [Biopython 1.87 dictionary](https://github.com/biopython/biopython/blob/biopython-187/Bio/Restriction/Restriction_Dictionary.py)
- [Biopython restriction documentation](https://biopython.org/DIST/docs/cookbook/Restriction.html)
- [NEB isoschizomers and neoschizomers](https://www.neb.com/en-ca/tools-and-resources/selection-charts/isoschizomers)
- [NEB Enzyme Finder](https://enzymefinder.neb.com/)
- [NEB methylation-dependent enzymes](https://www.neb.com/en-us/products/epigenetics/methylation-sensitive-restriction-enzymes/methylation-sensitive-restriction-enzymes)
