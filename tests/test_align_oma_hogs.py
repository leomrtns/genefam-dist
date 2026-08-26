import importlib.util
import signal
import shutil
import sys
from pathlib import Path

import pytest


SCRIPT = Path(__file__).parents[1] / "scripts" / "align_oma_hogs.py"
SPEC = importlib.util.spec_from_file_location("align_oma_hogs", SCRIPT)
aligner = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = aligner
SPEC.loader.exec_module(aligner)


AA_RECORDS = {
  "LACPS00918": "MSTQSPVHRNRVLTLIRSYYPNLSVTDRKIADYIIADPIKTAAQSISDLAAAVGVSTATV",
  "LACPS00919": "MSTQSPVHRNRVLTLIRSYYPNLSVTDRKIADYIIADPIKTAAQSISDLAAAVGVSTATA",
  "LACPS00920": "MSTQSPVHRNRVLTLIRSYYPNLSVTDRKIADYIIADPIKTAAQSISDLAAAVGVSTSAV",
  "LACPS00921": "MSTQSPVHRNRVLTLIRSYYPNLSVTDRKIADYIIADPIKTAAQSISDLAAAVGVSSATV",
}
THREE_DI_RECORDS = {
  "LACPS00918": "DPPPPPLPPDDLLVLLVVCLVPDDPLLNLLSVVCLVPLQVLLVDDLVVSCVVSVHDSSSD",
  "LACPS00919": "DPPPPPLPPDDLLVLLVVCLVPDDPLLNLLSVVCLVPLQVLLVDDLVVSCVVSVHDASSD",
  "LACPS00920": "DPPPPPLPPDDLLVLLVVCLVPDDPLLNLLSVVCLVPLQVLLVDDLVVSCVVSVHDSVSD",
  "LACPS00921": "DPPPPPLPPDDLLVLLVVCLVPDDPLLNLLSVVCLVPLQVLLVDDLVVSCVVSVHSSVSD",
}


def fasta_text(records):
  return "".join(f">{identifier}\n{sequence}\n" for identifier, sequence in records.items())


def write_family(tmp_path):
  family_dir = tmp_path / "families"
  member_dir = tmp_path / "members"
  family_dir.mkdir()
  member_dir.mkdir()
  (family_dir / "HOG_F0981989.faa").write_text(fasta_text(AA_RECORDS), encoding="utf-8")
  (family_dir / "HOG_F0981989.3di.fasta").write_text(
    fasta_text(THREE_DI_RECORDS),
    encoding="utf-8",
  )
  (member_dir / "HOG_F0981989.members.tsv").write_text(
    "omaid\tcanonical_id\nLACPS00918\tWP_013355364.1\n",
    encoding="utf-8",
  )
  return aligner.discover_families(tmp_path)[0]


def installed_foldmason():
  adjacent = Path(sys.executable).with_name("foldmason")
  return shutil.which("foldmason") or (str(adjacent) if adjacent.is_file() else None)


def test_family_inputs_accept_realistic_refseq_oma_members(tmp_path):
  family = write_family(tmp_path)

  amino_acids, three_di = aligner.validate_family(family)

  assert amino_acids == AA_RECORDS
  assert three_di == THREE_DI_RECORDS


def test_family_inputs_require_positionally_paired_sequences(tmp_path):
  family = write_family(tmp_path)
  mismatched = dict(THREE_DI_RECORDS)
  mismatched["LACPS00918"] = mismatched["LACPS00918"][:-1]
  family.three_di.write_text(fasta_text(mismatched), encoding="utf-8")

  with pytest.raises(aligner.AlignmentError, match="AA/3Di length mismatch"):
    aligner.validate_family(family)


def test_real_foldmason_aligns_oma_aa_and_3di_databases(tmp_path):
  foldmason = installed_foldmason()
  if foldmason is None:
    pytest.skip("FoldMason is not installed")
  write_family(tmp_path)
  output_dir = tmp_path / "output"

  result = aligner.main(
    [
      "--input-dir",
      str(tmp_path),
      "--output-dir",
      str(output_dir),
      "--foldmason",
      foldmason,
      "--threads",
      "1",
    ]
  )

  assert result == 0
  prefix = output_dir / "alignments" / "HOG_F0981989"
  aligned_aa = aligner.read_fasta(Path(f"{prefix}_aa.fa"))
  aligned_3di = aligner.read_fasta(Path(f"{prefix}_3di.fa"))
  assert aligned_aa.keys() == AA_RECORDS.keys()
  assert {len(sequence) for sequence in aligned_aa.values()} == {
    len(next(iter(aligned_3di.values())))
  }
  for omaid in AA_RECORDS:
    assert aligned_aa[omaid].replace("-", "") == AA_RECORDS[omaid]
    assert aligned_3di[omaid].replace("-", "") == THREE_DI_RECORDS[omaid]


def test_completed_alignment_is_a_checkpoint(tmp_path):
  foldmason = installed_foldmason()
  if foldmason is None:
    pytest.skip("FoldMason is not installed")
  write_family(tmp_path)
  output_dir = tmp_path / "output"
  arguments = [
    "--input-dir",
    str(tmp_path),
    "--output-dir",
    str(output_dir),
    "--foldmason",
    foldmason,
    "--threads",
    "1",
  ]

  assert aligner.main(arguments) == 0
  assert aligner.main(arguments) == 0

  manifest = (output_dir / "structure-alignments.tsv").read_text(encoding="utf-8")
  assert "HOG_F0981989\t4\texisting" in manifest


def test_foldmason_segfault_is_retried_with_one_thread(tmp_path, monkeypatch):
  family = write_family(tmp_path)
  amino_acids, three_di = aligner.validate_family(family)
  output_dir = tmp_path / "output"
  structure_commands = []

  def fake_run(command, started_at, label):
    if "structuremsa" not in command:
      return
    structure_commands.append(command.copy())
    if len(structure_commands) == 1:
      raise aligner.CommandError(command, -signal.SIGSEGV)
    prefix = Path(command[3])
    Path(f"{prefix}_aa.fa").write_text(fasta_text(amino_acids), encoding="utf-8")
    Path(f"{prefix}_3di.fa").write_text(fasta_text(three_di), encoding="utf-8")
    Path(f"{prefix}.nw").write_text("(LACPS00918,LACPS00919);\n", encoding="utf-8")

  monkeypatch.setattr(aligner, "run_command", fake_run)
  aligner.align_family(
    "foldmason",
    family,
    amino_acids,
    three_di,
    output_dir,
    8,
    0,
    False,
    0,
  )

  assert structure_commands[0][structure_commands[0].index("--threads") + 1] == "8"
  assert structure_commands[1][structure_commands[1].index("--threads") + 1] == "1"


def test_failure_details_are_written_only_to_error_report(tmp_path, monkeypatch):
  write_family(tmp_path)
  output_dir = tmp_path / "output"
  detail = "Command failed with exit status -11: a very long command"
  monkeypatch.setattr(aligner, "find_foldmason", lambda command: command)
  monkeypatch.setattr(
    aligner,
    "align_family",
    lambda *args, **kwargs: (_ for _ in ()).throw(aligner.AlignmentError(detail)),
  )

  assert aligner.main(["--input-dir", str(tmp_path), "--output-dir", str(output_dir)]) == 1

  manifest = (output_dir / "structure-alignments.tsv").read_text(encoding="utf-8")
  errors = (output_dir / "structure-alignment-errors.tsv").read_text(encoding="utf-8")
  assert "HOG_F0981989\t4\terror" in manifest
  assert detail not in manifest
  assert detail in errors
