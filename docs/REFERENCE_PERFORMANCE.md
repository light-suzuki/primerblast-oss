# Read-only reference access

GFF gene lookup and locus annotations now use an on-demand SQLite query cache.
The first request streams the GFF (including gzip) once and constructs a derived
index. Subsequent requests load only the selected genes and their descendants,
preserving aliases, ambiguity checks, feature ordering, and full transcript models.
An interval with no genes remains different from an unknown chromosome.

The source GFF, FASTA, `.fai`, and BLAST databases are opened for reading only.
No indexes are created alongside references. The disposable query cache is in
the user's `.cache/sequence-reference-queries` directory, or the separate folder
configured by `SEQWB_REFERENCE_CACHE_DIR`. It contains annotation data and must
remain private, like the source reference. Do not upload it. Cache fingerprints
include the source path, size and modification/change timestamps; changed inputs
produce a new cache. Only completed builds are published atomically. Old cache
files may be removed while the app is stopped; they will be regenerated.

FASTA access continues to use existing `.fai` files, reads exactly the requested
byte range, and reuses up to eight small index dictionaries with timestamp-based
invalidation. The search engine, mismatch rules and scientific cutoffs are unchanged.

Run `python scripts/benchmark_reference_access.py` for a bounded synthetic
benchmark with 5,000 genes / 15,000 features. It compares whole-file parsing,
initial index construction and a warm lookup, checks identical exon coordinates,
and verifies the source hash did not change. Timing and Python allocation peaks
are machine-specific; they are not a real-genome BLAST benchmark or Wet evidence.
The initial build can be slower and uses additional cache disk space. This change
does not claim to accelerate the BLAST search algorithm itself.
