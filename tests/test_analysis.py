import csv
import math
import subprocess
import sys

import pytest

from treesignal import TreeSpace
from treesignal.analysis import TreeAnalysis, analyze_files, generate_references


FIRST = "((A:1,B:1):1,(C:1,D:1):1);"
SECOND = "((A:1,C:1):1,(B:1,D:1):1);"
SCALED = "((A:2,B:2):2,(C:2,D:2):2);"


def test_length_metrics_rank_the_matching_reference_first():
    references = TreeSpace(FIRST + SECOND, rooted=True)
    analysis = TreeAnalysis(references)

    result = analysis.compare(FIRST)

    assert list(result) == ["patristic", "wrf", "kf"]
    for values in result.values():
        assert values[0] == 0
        assert values[1] > 0
    assert analysis.references is references


def test_scale_normalization_removes_uniform_rate_changes():
    analysis = TreeAnalysis(FIRST)

    raw = analysis.compare(SCALED)
    normalized = analysis.compare(SCALED, normalize=True)

    assert all(values[0] > 0 for values in raw.values())
    assert all(values[0] == pytest.approx(0) for values in normalized.values())


def test_reinitialized_references_cannot_leave_stale_matrix_results():
    references = TreeSpace(FIRST)
    analysis = TreeAnalysis(references)
    assert analysis.compare(FIRST, metrics="patristic")["patristic"] == (0,)

    references.__init__(SCALED)

    assert analysis.compare(FIRST, metrics="patristic")["patristic"][0] > 0


def test_paralog_patristic_projection_uses_explicit_mapping():
    gene = "((A1:1,A2:3):1,(B1:1,C1:1):1);"
    reference = "(A:1,(B:1,C:1):1);"
    mapping = {"A1": "A", "A2": "A", "B1": "B", "C1": "C"}
    analysis = TreeAnalysis(reference)

    minimum = analysis.compare(gene, metrics="patristic", gene_to_species=mapping)
    mean = analysis.compare(gene, metrics="patristic", gene_to_species=mapping, aggregation="mean")

    assert minimum["patristic"] == pytest.approx((math.sqrt(2 / 3),))
    assert mean["patristic"] == pytest.approx((math.sqrt(8 / 3),))
    with pytest.raises(ValueError, match="supported for patristic"):
        analysis.compare(gene, gene_to_species=mapping)


def test_patristic_restricts_to_shared_species():
    analysis = TreeAnalysis(FIRST)
    gene = "(A:1,(B:1,X:1):1);"

    assert analysis.compare(gene, metrics="patristic")["patristic"] == pytest.approx((1,))


@pytest.mark.parametrize("kwargs, message", [
    ({"metrics": []}, "at least one metric"),
    ({"metrics": ["unknown"]}, "unknown metrics"),
    ({"metrics": ["patristic", "patristic"]}, "duplicates"),
    ({"aggregation": "median"}, "aggregation"),
    ({"gene_to_species": {"A": "A"}}, "missing gene-to-species"),
])
def test_invalid_analysis_options_are_rejected(kwargs, message):
    with pytest.raises(ValueError, match=message):
        TreeAnalysis(FIRST).compare(FIRST, **kwargs)


def test_gene_file_must_contain_exactly_one_tree():
    with pytest.raises(ValueError, match="exactly one"):
        TreeAnalysis(FIRST).compare(FIRST + SECOND)


def test_missing_lengths_cannot_silently_become_unit_lengths():
    with pytest.raises(ValueError, match="require complete"):
        TreeAnalysis(FIRST).compare("((A,B),(C,D));", metrics="patristic")


def test_generate_references_is_seeded_and_uses_requested_names():
    labels = ("species_A", "species_B", "species_C", "species_D")
    first = generate_references(labels, count=5, seed=12)
    second = generate_references(labels, count=5, seed=12)

    assert first.newick == second.newick
    assert set(first.taxa) == set(labels)
    assert first.ntrees == first.ndistinct == 5
    assert first.has_branch_lengths


def test_generated_branch_scale_is_explicit():
    first = generate_references(("A", "B", "C", "D"), count=1, seed=12)
    second = generate_references(("A", "B", "C", "D"), count=1, seed=12, branch_length=2)

    matrix1, matrix2 = first.patristic_matrix(), second.patristic_matrix()
    assert matrix2[0][1] == pytest.approx(2 * matrix1[0][1])


@pytest.mark.parametrize("kwargs", [
    {"count": 0}, {"n_spr": 0}, {"seed": -1},
    {"branch_length": float("nan")}, {"branch_length": -1},
])
def test_invalid_reference_design_is_rejected(kwargs):
    with pytest.raises(ValueError):
        generate_references(("A", "B", "C", "D"), **kwargs)


def test_batch_output_has_explicit_metric_reference_and_family(tmp_path):
    (tmp_path / "b.nwk").write_text(SECOND)
    (tmp_path / "a.nwk").write_text(FIRST)
    output = tmp_path / "result.tsv"

    rows = analyze_files(str(tmp_path / "*.nwk"), FIRST + SECOND, output=output)

    assert len(rows) == 12
    assert rows[0] == {"family": "a.nwk", "reference": 1, "metric": "patristic", "value": 0}
    with output.open() as handle:
        saved = list(csv.DictReader(handle, delimiter="\t"))
    assert len(saved) == len(rows)
    assert saved[0]["family"] == "a.nwk"


def test_batch_failure_keeps_existing_output_and_identifies_family(tmp_path):
    (tmp_path / "bad.nwk").write_text("(A:1,B:1);")
    output = tmp_path / "result.tsv"
    output.write_text("previous results\n")

    with pytest.raises(ValueError, match="bad.nwk"):
        analyze_files([tmp_path / "bad.nwk"], FIRST, output=output)

    assert output.read_text() == "previous results\n"


def test_batch_refuses_to_overwrite_an_input_tree(tmp_path):
    source = tmp_path / "gene.nwk"
    source.write_text(FIRST)

    with pytest.raises(ValueError, match="overwrite"):
        analyze_files([source], FIRST, output=source)

    assert source.read_text() == FIRST


def test_cli_generates_references_and_analyzes_family(tmp_path):
    species = tmp_path / "species.txt"
    species.write_text("A\nB\nC\nD\n")
    refs = tmp_path / "refs.nwk"
    subprocess.run([sys.executable, "-m", "treesignal", "references", "--species", str(species),
                    "--count", "3", "--seed", "10", "--output", str(refs)], check=True)
    genes = tmp_path / "genes"
    genes.mkdir()
    (genes / "family.nwk").write_text(FIRST)
    output = tmp_path / "result.tsv"
    subprocess.run([sys.executable, "-m", "treesignal", "analyze", "--references", str(refs),
                    "--genes", str(genes / "*.nwk"), "--output", str(output)], check=True)

    with output.open() as handle:
        assert len(list(csv.DictReader(handle, delimiter="\t"))) == 9


def test_cli_mapping_file_supports_paralogs(tmp_path):
    (tmp_path / "refs.nwk").write_text("(A:1,(B:1,C:1):1);")
    (tmp_path / "gene.nwk").write_text("((A1:1,A2:3):1,(B1:1,C1:1):1);")
    (tmp_path / "mapping.tsv").write_text("A1\tA\nA2\tA\nB1\tB\nC1\tC\n")
    result = subprocess.run([sys.executable, "-m", "treesignal", "analyze",
                             "--references", str(tmp_path / "refs.nwk"),
                             "--genes", str(tmp_path / "gene.nwk"),
                             "--mapping", str(tmp_path / "mapping.tsv"),
                             "--metrics", "patristic", "--output", str(tmp_path / "result.tsv")],
                            capture_output=True, text=True)

    assert result.returncode == 0, result.stderr
