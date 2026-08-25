import importlib.util
import json
import sys
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "scripts" / "align_oma_hogs.py"
SPEC = importlib.util.spec_from_file_location("align_oma_hogs", SCRIPT)
aligner = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = aligner
SPEC.loader.exec_module(aligner)


def write_family(tmp_path):
  family_dir = tmp_path / "families"
  member_dir = tmp_path / "members"
  family_dir.mkdir()
  member_dir.mkdir()
  (family_dir / "HOG_F1.faa").write_text(
    ">OMA1\nMPEPTIDE\n>OMA2\nMSEQUENCE\n",
    encoding="utf-8",
  )
  (family_dir / "HOG_F1.3di.fasta").write_text(
    ">OMA1\nACDEFGHI\n>OMA2\nKLMNPQRSQ\n",
    encoding="utf-8",
  )
  (member_dir / "HOG_F1.members.tsv").write_text(
    "omaid\tcanonical_id\nOMA1\tP00001\nOMA2\tP00002\n",
    encoding="utf-8",
  )
  return aligner.discover_families(tmp_path)[0]


def test_family_inputs_have_parallel_amino_acid_and_3di_sequences(tmp_path):
  family = write_family(tmp_path)

  amino_acids, three_di = aligner.validate_family_sequences(family)

  assert amino_acids == {"OMA1": "MPEPTIDE", "OMA2": "MSEQUENCE"}
  assert three_di == {"OMA1": "ACDEFGHI", "OMA2": "KLMNPQRSQ"}


def test_select_model_requires_exact_oma_sequence():
  records = [
    {"sequence": "OLD", "cifUrl": "https://example.test/old.cif", "latestVersion": 5},
    {"sequence": "MPEPTIDE", "cifUrl": "https://example.test/new.cif", "latestVersion": 6},
  ]

  selected = aligner.select_model(records, "P00001", "MPEPTIDE")

  assert selected["cifUrl"] == "https://example.test/new.cif"


def test_cached_oma_json_supplies_missing_canonical_id(tmp_path):
  family = write_family(tmp_path)
  assert family.members is not None
  family.members.unlink()
  family = aligner.discover_families(tmp_path)[0]
  for omaid, canonical, sequence in (
    ("OMA1", "P00001", "MPEPTIDE"),
    ("OMA2", "P00002", "MSEQUENCE"),
  ):
    cache_path = tmp_path / ".cache" / omaid / "record.json"
    cache_path.parent.mkdir(parents=True)
    cache_path.write_text(
      json.dumps({"omaid": omaid, "canonicalid": canonical, "sequence": sequence}),
      encoding="utf-8",
    )
  amino_acids, _three_di = aligner.validate_family_sequences(family)

  result = aligner.canonical_ids(
    family,
    amino_acids,
    {},
    aligner.cached_oma_proteins(tmp_path / ".cache"),
  )

  assert result == {"OMA1": "P00001", "OMA2": "P00002"}


def test_validate_foldmason_output_checks_projected_amino_acids(tmp_path):
  family = write_family(tmp_path)
  output_prefix = tmp_path / "alignments" / family.stem
  output_prefix.parent.mkdir()
  Path(f"{output_prefix}_aa.fa").write_text(
    ">OMA1\nMP--EPTIDE\n>OMA2\nM-SEQUENCE\n",
    encoding="utf-8",
  )
  Path(f"{output_prefix}_3di.fa").write_text(
    ">OMA1\nAC--DEFGHI\n>OMA2\nK-LMNPQRSQ\n",
    encoding="utf-8",
  )

  aligner.validate_foldmason_output(
    family,
    output_prefix,
    {"OMA1": "MPEPTIDE", "OMA2": "MSEQUENCE"},
  )


def test_download_family_models_is_restartable(tmp_path):
  family = write_family(tmp_path)
  amino_acids, _three_di = aligner.validate_family_sequences(family)

  class FakeClient:
    refresh = False

    def __init__(self):
      self.model_calls = 0

    def metadata(self, accession):
      sequence = amino_acids["OMA1" if accession == "P00001" else "OMA2"]
      return [{"sequence": sequence, "cifUrl": f"https://example.test/{accession}.cif"}]

    def model(self, url, path):
      if path.is_file() and path.stat().st_size > 0:
        return
      self.model_calls += 1
      aligner.atomic_write_text(path, "data_model\n")

  client = FakeClient()
  identifiers = {"OMA1": "P00001", "OMA2": "P00002"}
  structure_dir = tmp_path / "structures"

  first = aligner.download_family_models(
    client,
    family,
    amino_acids,
    identifiers,
    structure_dir,
    workers=2,
    started_at=aligner.time.monotonic(),
  )
  second = aligner.download_family_models(
    client,
    family,
    amino_acids,
    identifiers,
    structure_dir,
    workers=2,
    started_at=aligner.time.monotonic(),
  )

  assert first == second
  assert client.model_calls == 2


def test_run_foldmason_stages_and_publishes_outputs(tmp_path):
  family = write_family(tmp_path)
  structure_dir = tmp_path / "structures"
  structure_dir.mkdir()
  output_prefix = tmp_path / "output" / "alignments" / family.stem
  original_run = aligner.subprocess.run

  def fake_run(command, check):
    assert check
    staged_prefix = Path(command[3])
    Path(f"{staged_prefix}_aa.fa").write_text(
      ">OMA1\nMP--EPTIDE\n>OMA2\nM-SEQUENCE\n",
      encoding="utf-8",
    )
    Path(f"{staged_prefix}_3di.fa").write_text(
      ">OMA1\nAC--DEFGHI\n>OMA2\nK-LMNPQRSQ\n",
      encoding="utf-8",
    )
    Path(f"{staged_prefix}.nw").write_text("(OMA1,OMA2);\n", encoding="utf-8")

  aligner.subprocess.run = fake_run
  try:
    aligner.run_foldmason(
      "foldmason",
      family,
      structure_dir,
      output_prefix,
      {"OMA1": "MPEPTIDE", "OMA2": "MSEQUENCE"},
      report_mode=0,
    )
  finally:
    aligner.subprocess.run = original_run

  assert Path(f"{output_prefix}_aa.fa").is_file()
  assert Path(f"{output_prefix}_3di.fa").is_file()
  assert Path(f"{output_prefix}.nw").is_file()


def test_download_only_main_uses_downloader_output(tmp_path):
  write_family(tmp_path)
  output_dir = tmp_path / "structure-output"
  sequences = {"P00001": "MPEPTIDE", "P00002": "MSEQUENCE"}
  original_client = aligner.AFDBClient

  class FakeClient:
    refresh = False

    def __init__(self, *args):
      pass

    def metadata(self, accession):
      return [
        {
          "sequence": sequences[accession],
          "cifUrl": f"https://example.test/{accession}.cif",
        }
      ]

    def model(self, url, path):
      aligner.atomic_write_text(path, "data_model\n")

  aligner.AFDBClient = FakeClient
  try:
    result = aligner.main(
      [
        "--input-dir",
        str(tmp_path),
        "--output-dir",
        str(output_dir),
        "--download-only",
      ]
    )
  finally:
    aligner.AFDBClient = original_client

  assert result == 0
  manifest = (output_dir / "structure-alignments.tsv").read_text(encoding="utf-8")
  assert "HOG_F1\t2\tdownloaded" in manifest
  assert len(list((output_dir / "structures" / "HOG_F1").glob("*.cif"))) == 2
