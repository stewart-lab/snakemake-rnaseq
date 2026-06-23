"""Unit tests for app/workflow/scripts/gtf_gene_map.py."""
import csv
import gzip

import pytest

import gtf_gene_map


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def gtf_line(attributes, feature="gene", seqname="1"):
    """A full 9-column GTF/GFF line with the given attributes column."""
    cols = [seqname, "src", feature, "1", "100", ".", "+", ".", attributes]
    return "\t".join(cols)


def write_lines(tmp_path, lines, name="ann.gtf"):
    path = tmp_path / name
    path.write_text("\n".join(lines) + "\n")
    return str(path)


# --------------------------------------------------------------------------- #
# build_map: per-source gene name extraction
# --------------------------------------------------------------------------- #
def test_ensembl_uses_gene_name(tmp_path):
    path = write_lines(tmp_path, [
        gtf_line('gene_id "ENSG001"; gene_version "5"; gene_name "TP53"; gene_biotype "protein_coding";'),
    ])
    assert gtf_gene_map.build_map(path) == {"ENSG001": "TP53"}


def test_ncbi_uses_gene_attribute(tmp_path):
    path = write_lines(tmp_path, [
        gtf_line('gene_id "DDX11L1"; transcript_id ""; gene "DDX11L1"; gene_biotype "misc_RNA";'),
    ])
    assert gtf_gene_map.build_map(path) == {"DDX11L1": "DDX11L1"}


def test_gene_symbol_attribute(tmp_path):
    path = write_lines(tmp_path, [gtf_line('gene_id "G42"; gene_symbol "FOO1";')])
    assert gtf_gene_map.build_map(path) == {"G42": "FOO1"}


def test_gff3_name_attribute(tmp_path):
    path = write_lines(tmp_path, [
        gtf_line("ID=gene:ENSG002;Name=BRCA1;biotype=protein_coding;gene_id=ENSG002"),
    ])
    assert gtf_gene_map.build_map(path) == {"ENSG002": "BRCA1"}


# --------------------------------------------------------------------------- #
# build_map: tricky cases
# --------------------------------------------------------------------------- #
def test_havana_gene_is_not_mistaken_for_a_name(tmp_path):
    # "havana_gene" must NOT satisfy the "gene" key (exact-match only); with no
    # real name attribute the gene falls back to its id.
    path = write_lines(tmp_path, [
        gtf_line('gene_id "ENSG009"; havana_gene "OTTHUMG00000001";'),
    ])
    assert gtf_gene_map.build_map(path) == {"ENSG009": "ENSG009"}


def test_gene_biotype_is_not_mistaken_for_a_name(tmp_path):
    path = write_lines(tmp_path, [gtf_line('gene_id "G1"; gene_biotype "lincRNA";')])
    assert gtf_gene_map.build_map(path) == {"G1": "G1"}


def test_nameless_gene_falls_back_to_gene_id(tmp_path):
    path = write_lines(tmp_path, [gtf_line('gene_id "G99";')])
    assert gtf_gene_map.build_map(path) == {"G99": "G99"}


def test_gene_name_takes_priority_over_gene(tmp_path):
    path = write_lines(tmp_path, [
        gtf_line('gene_id "G1"; gene "OTHER"; gene_name "REAL";'),
    ])
    assert gtf_gene_map.build_map(path) == {"G1": "REAL"}


def test_fallback_is_upgraded_to_a_real_name_on_a_later_line(tmp_path):
    # First line has no name (fallback), a later line for the same gene does.
    path = write_lines(tmp_path, [
        gtf_line('gene_id "G1";', feature="gene"),
        gtf_line('gene_id "G1"; gene_name "SYM";', feature="transcript"),
    ])
    assert gtf_gene_map.build_map(path) == {"G1": "SYM"}


def test_existing_name_is_not_downgraded_by_a_later_nameless_line(tmp_path):
    path = write_lines(tmp_path, [
        gtf_line('gene_id "G1"; gene_name "SYM";', feature="gene"),
        gtf_line('gene_id "G1";', feature="exon"),
    ])
    assert gtf_gene_map.build_map(path) == {"G1": "SYM"}


def test_gene_id_version_suffix_is_preserved(tmp_path):
    path = write_lines(tmp_path, [gtf_line('gene_id "ENSG002.3"; gene_name "BRCA1";')])
    assert gtf_gene_map.build_map(path) == {"ENSG002.3": "BRCA1"}


def test_comment_and_short_lines_are_skipped(tmp_path):
    path = write_lines(tmp_path, [
        "#!genome-build GRCh38",
        "## a comment",
        "1\tsrc\tgene\t1\t100",  # too few columns
        gtf_line('gene_id "G1"; gene_name "TP53";'),
    ])
    assert gtf_gene_map.build_map(path) == {"G1": "TP53"}


def test_insertion_order_is_preserved(tmp_path):
    path = write_lines(tmp_path, [
        gtf_line('gene_id "C"; gene_name "c";'),
        gtf_line('gene_id "A"; gene_name "a";'),
        gtf_line('gene_id "B"; gene_name "b";'),
    ])
    assert list(gtf_gene_map.build_map(path).keys()) == ["C", "A", "B"]


def test_gzipped_input_is_read(tmp_path):
    path = tmp_path / "ann.gtf.gz"
    with gzip.open(path, "wt") as fh:
        fh.write(gtf_line('gene_id "G1"; gene_name "TP53";') + "\n")
    assert gtf_gene_map.build_map(str(path)) == {"G1": "TP53"}


def test_mixed_sources_in_one_file(tmp_path):
    path = write_lines(tmp_path, [
        gtf_line('gene_id "ENSG001"; gene_name "TP53";'),
        gtf_line('gene_id "DDX11L1"; gene "DDX11L1";'),
        gtf_line('gene_id "G42"; gene_symbol "FOO1";'),
        gtf_line("gene_id=G7;Name=WIDGET1"),
        gtf_line('gene_id "G99"; gene_biotype "lincRNA";'),
    ])
    assert gtf_gene_map.build_map(path) == {
        "ENSG001": "TP53",
        "DDX11L1": "DDX11L1",
        "G42": "FOO1",
        "G7": "WIDGET1",
        "G99": "G99",
    }


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def test_main_writes_expected_csv(tmp_path, monkeypatch):
    gtf = write_lines(tmp_path, [
        gtf_line('gene_id "ENSG001"; gene_name "TP53";'),
        gtf_line('gene_id "G99";'),
    ])
    out = tmp_path / "gene_map.csv"
    monkeypatch.setattr("sys.argv", ["gtf_gene_map.py", "--gtf", gtf, "--out", str(out)])

    gtf_gene_map.main()

    with open(out, newline="") as fh:
        rows = list(csv.reader(fh))
    assert rows[0] == ["gene_id", "gene_name"]
    assert rows[1:] == [["ENSG001", "TP53"], ["G99", "G99"]]


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
