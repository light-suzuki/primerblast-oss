"""Tests for versioned and optional strong provenance fingerprints."""
import hashlib

from primerblast_oss import provenance


def test_model_versions_are_explicit_and_independent():
    versions = provenance.model_versions()
    assert provenance.PROVENANCE_SCHEMA_VERSION
    assert set(versions) == {"specificity", "gel", "risk", "restriction"}
    assert all(isinstance(value, str) and value for value in versions.values())


def test_lightweight_fingerprint_does_not_compute_full_sha256(tmp_path, monkeypatch):
    path = tmp_path / "reference.fa"
    path.write_text(">chr1\nACGTACGT\n", encoding="utf-8")

    def forbidden(*_args, **_kwargs):
        raise AssertionError("lightweight provenance must not full-hash the file")

    monkeypatch.setattr(provenance, "sha256_file", forbidden)
    result = provenance.fasta_fingerprint(str(path), strong=False)
    assert result["sha1_prefix"]
    assert result["sha256"] is None


def test_strong_file_fingerprint_is_full_sha256(tmp_path):
    path = tmp_path / "reference.fa"
    payload = b">chr1\nACGTACGT\n"
    path.write_bytes(payload)
    result = provenance.fasta_fingerprint(str(path), strong=True)
    assert result["sha256"] == hashlib.sha256(payload).hexdigest()


def test_blast_database_strong_hashes_all_nucleotide_index_files(tmp_path):
    prefix = tmp_path / "db"
    files = {
        ".nin": b"index-in",
        ".nsq": b"index-seq",
        ".nhr": b"index-header",
        ".fai": b"not-a-blast-index",
    }
    for suffix, payload in files.items():
        (tmp_path / ("db" + suffix)).write_bytes(payload)

    result = provenance.db_fingerprint(str(prefix), strong=True)
    names = {entry["name"] for entry in result["files"]}
    assert names == {"db.nhr", "db.nin", "db.nsq"}
    by_name = {entry["name"]: entry for entry in result["files"]}
    assert by_name["db.nin"]["sha256"] == hashlib.sha256(b"index-in").hexdigest()


def test_manifest_records_schema_models_revision_and_reference_policy(
        tmp_path, monkeypatch):
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
    db = tmp_path / "genome_db"
    (tmp_path / "genome_db.nin").write_bytes(b"nin")
    (tmp_path / "genome_db.nsq").write_bytes(b"nsq")
    (tmp_path / "genome_db.nhr").write_bytes(b"nhr")
    monkeypatch.setenv("PRIMERBLAST_GIT_REVISION", "abc123")
    monkeypatch.setattr(
        provenance, "tool_versions", lambda: {"primerblast_oss": "test"})

    manifest = provenance.make_manifest(
        {"example": True}, [str(db)],
        thermo_genomes={str(db): str(fasta)},
        reference_files={"design_genome": str(fasta)},
        strong_hashes=True,
    )
    assert manifest["provenance_schema_version"] == provenance.PROVENANCE_SCHEMA_VERSION
    assert manifest["source_revision"] == "abc123"
    assert manifest["model_versions"] == provenance.model_versions()
    assert manifest["fingerprint_policy"]["mode"] == "full_sha256"
    assert manifest["thermo_genomes"][str(db)]["sha256"]
    assert manifest["reference_files"]["design_genome"]["sha256"]
    assert all(entry["sha256"] for entry in manifest["databases"][0]["files"])
