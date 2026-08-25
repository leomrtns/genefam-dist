import csv
import importlib.util
import sys
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "scripts" / "download_oma_hogs.py"
SPEC = importlib.util.spec_from_file_location("download_oma_hogs", SCRIPT)
oma = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = oma
SPEC.loader.exec_module(oma)


class FakeClient:
    def __init__(self, responses):
        self.responses = responses
        self.calls = []

    def get_json(self, path, params=None):
        self.calls.append((path, params))
        key = (path, tuple(sorted((params or {}).items())))
        return self.responses[key]


def test_iter_hogs_at_level_uses_all_pages():
    client = FakeClient(
        {
            ("/hog/", (("level", "Fungi"), ("page", 1), ("per_page", 2))): [
                {"hog_id": "HOG:F0000001"},
                {"hog_id": "HOG:F0000002"},
            ],
            ("/hog/", (("level", "Fungi"), ("page", 2), ("per_page", 2))): [
                {"hog_id": "HOG:F0000003"}
            ],
        }
    )

    assert [hog["hog_id"] for hog in oma.iter_hogs_at_level(client, "Fungi", 2)] == [
        "HOG:F0000001",
        "HOG:F0000002",
        "HOG:F0000003",
    ]


def test_root_grouping_keeps_ancestral_gene_components():
    hogs = [
        {"hog_id": "HOG:F1029536.1a", "roothog_id": 1029536},
        {"hog_id": "HOG:F1029536.1b", "roothog_id": 1029536},
        {"hog_id": "HOG:F1029537", "roothog_id": 1029537},
    ]

    families = list(oma.group_hogs(hogs, "root-hog"))

    assert [family.family_id for family in families] == [
        "HOG:F1029536",
        "HOG:F1029537",
    ]
    assert [component["hog_id"] for component in families[0].components] == [
        "HOG:F1029536.1a",
        "HOG:F1029536.1b",
    ]


def test_root_filter_keeps_low_completeness_components_after_family_selected():
    family = oma.Family(
        "HOG:F1029536",
        "HOG:F1029536",
        (
            {"hog_id": "HOG:F1029536.1a", "completeness_score": 0.8, "nr_genes": 4},
            {"hog_id": "HOG:F1029536.1b", "completeness_score": 0.1, "nr_genes": 2},
        ),
    )

    assert oma.family_metadata_filter(family, "root-hog", 0.3, 4, None)
    assert not oma.family_metadata_filter(family, "ancestral-gene", 0.9, 4, None)


def test_download_merged_family_records_component_assignment(tmp_path):
    components = (
        {
            "hog_id": "HOG:F1029536.1a",
            "roothog_id": 1029536,
            "completeness_score": 0.8,
            "description": "example family",
        },
        {
            "hog_id": "HOG:F1029536.1b",
            "roothog_id": 1029536,
            "completeness_score": 0.2,
            "description": "example family",
        },
    )
    members = {}
    proteins = {}
    for index, (hog_id, code) in enumerate(
        [
            ("HOG:F1029536.1a", "SP001"),
            ("HOG:F1029536.1a", "SP002"),
            ("HOG:F1029536.1b", "SP003"),
            ("HOG:F1029536.1b", "SP004"),
        ],
        start=1,
    ):
        member = {
            "entry_nr": index,
            "omaid": f"{code}00001",
            "oma_hog_id": hog_id,
            "sequence_md5": f"md5-{index}",
            "species": {"code": code, "species": f"Species {index}", "taxon_id": index},
        }
        members.setdefault(hog_id, []).append(member)
        proteins[index] = {
            "entry_nr": index,
            "omaid": member["omaid"],
            "canonicalid": f"CAN{index}",
            "sequence": "MPEPTIDE",
            "cdna": "ATGCCCGAACTGACCATCGATGAA",
            "structure": {"source": "AlphaFold", "sequence_3di": "ACDEFGHI"},
            "xref": f"https://omabrowser.org/api/protein/{index}/xref/",
            "sequence_md5": f"md5-{index}",
        }

    responses = {}
    for hog_id, records in members.items():
        path = f"/hog/{hog_id}/members/"
        responses[(path, (("level", "o__Enterobacterales"),))] = {
            "hog_id": hog_id,
            "level": "o__Enterobacterales",
            "members": records,
        }
    for entry_nr, protein in proteins.items():
        responses[(f"/protein/{entry_nr}/", ())] = protein
        responses[(f"/protein/{entry_nr}/xref/", ())] = [
            {
                "xref": f"P0000{entry_nr}",
                "source": "UniProtKB/SwissProt",
                "seq_match": "exact",
            },
            {
                "xref": f"PROT{entry_nr}_SPECIES",
                "source": "UniProtKB/SwissProt",
                "seq_match": "exact",
            },
            {
                "xref": f"WP_00000000{entry_nr}.1",
                "source": "RefSeq",
                "seq_match": "exact",
            },
        ]
    client = FakeClient(responses)
    family = oma.Family("HOG:F1029536", "HOG:F1029536", components)

    row = oma.download_family(
        client,
        family,
        "o__Enterobacterales",
        tmp_path,
        "all",
        workers=2,
        min_members=4,
        min_species=4,
        max_members=None,
        force=False,
        include_xrefs=True,
    )

    assert row is not None
    fasta = (tmp_path / row["fasta"]).read_text()
    assert fasta.count(">") == 4
    assert (tmp_path / row["protein_fasta"]).read_text().count("MPEPTIDE") == 4
    assert (tmp_path / row["cdna_fasta"]).read_text().count("ATGCCCGAACTGACCATCGATGAA") == 4
    assert (tmp_path / row["three_di_fasta"]).read_text().count("ACDEFGHI") == 4
    with (tmp_path / row["members_tsv"]).open(newline="") as handle:
        metadata = list(csv.DictReader(handle, delimiter="\t"))
    assert {record["ancestor_hog_id"] for record in metadata} == {
        "HOG:F1029536.1a",
        "HOG:F1029536.1b",
    }
    assert {record["root_hog_id"] for record in metadata} == {"HOG:F1029536"}
    assert {record["uniprot_accessions"] for record in metadata} == {
        "P00001",
        "P00002",
        "P00003",
        "P00004",
    }
    assert {record["structure_source"] for record in metadata} == {"AlphaFold"}
    assert {record["xrefs_queried"] for record in metadata} == {"yes"}


def test_expected_member_count_sums_components():
    family = oma.Family(
        "HOG:F1029536",
        "HOG:F1029536",
        ({"hog_id": "HOG:F1029536.1a", "nr_genes": 4},
         {"hog_id": "HOG:F1029536.1b", "nr_genes": 2}),
    )

    assert oma.expected_member_count(family) == 6
