#!/usr/bin/env python3
"""Align OMA gene families from paired amino-acid and ProstT5 3Di FASTA files.

The input is produced by ``download_oma_hogs.py --sequence-type all``. For
each family, this script builds matching Foldseek amino-acid and 3Di databases
and runs FoldMason's ``structuremsa`` command. No coordinates, AlphaFold DB
accessions, network access, or structure prediction are required.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import signal
import shutil
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence


class AlignmentError(RuntimeError):
  """An input or FoldMason result was invalid."""


class CommandError(AlignmentError):
  """An external command returned a non-zero exit status."""

  def __init__(self, command: list[str], return_code: int):
    self.command = command
    self.return_code = return_code
    super().__init__(
      f"Command failed with exit status {return_code}: {' '.join(command)}"
    )


@dataclass(frozen=True)
class FamilyInput:
  stem: str
  amino_acids: Path
  three_di: Path


MANIFEST_FIELDS = (
  "family",
  "member_count",
  "status",
  "amino_acid_alignment",
  "three_di_alignment",
  "guide_tree",
)

ERROR_FIELDS = ("family", "member_count", "error")


def elapsed_minutes(started_at: float) -> float:
  return (time.monotonic() - started_at) / 60


def progress(started_at: float, message: str) -> None:
  print(
    f"[{elapsed_minutes(started_at):.1f} min] {message}",
    file=sys.stderr,
    flush=True,
  )


def atomic_write(path: Path, data: bytes) -> None:
  path.parent.mkdir(parents=True, exist_ok=True)
  fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
  try:
    with os.fdopen(fd, "wb") as handle:
      handle.write(data)
    os.replace(temporary, path)
  except BaseException:
    try:
      os.unlink(temporary)
    except FileNotFoundError:
      pass
    raise


def atomic_write_text(path: Path, text: str) -> None:
  atomic_write(path, text.encode("utf-8"))


def safe_name(value: str) -> str:
  return re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("._")


def read_fasta(path: Path) -> dict[str, str]:
  records: dict[str, list[str]] = {}
  identifier: str | None = None
  try:
    with path.open(encoding="utf-8") as handle:
      for line_number, raw_line in enumerate(handle, start=1):
        line = raw_line.strip()
        if not line:
          continue
        if line.startswith(">"):
          identifier = line[1:].split(None, 1)[0]
          if not identifier:
            raise AlignmentError(f"Empty FASTA identifier in {path}:{line_number}")
          if identifier in records:
            raise AlignmentError(f"Duplicate FASTA identifier {identifier} in {path}")
          records[identifier] = []
        elif identifier is None:
          raise AlignmentError(f"Sequence before first FASTA header in {path}:{line_number}")
        else:
          records[identifier].append(line)
  except OSError as error:
    raise AlignmentError(f"Could not read {path}: {error}") from error
  if not records:
    raise AlignmentError(f"No FASTA records in {path}")
  return {key: "".join(parts) for key, parts in records.items()}


def discover_families(input_dir: Path) -> list[FamilyInput]:
  family_dir = input_dir / "families"
  families: list[FamilyInput] = []
  for amino_acids in sorted(family_dir.glob("*.faa")):
    stem = amino_acids.name.removesuffix(".faa")
    three_di = family_dir / f"{stem}.3di.fasta"
    if not three_di.is_file():
      raise AlignmentError(f"Missing companion 3Di FASTA for {amino_acids}: {three_di}")
    families.append(FamilyInput(stem, amino_acids, three_di))
  if not families:
    raise AlignmentError(f"No family .faa files found under {family_dir}")
  return families


def validate_family(family: FamilyInput) -> tuple[dict[str, str], dict[str, str]]:
  amino_acids = read_fasta(family.amino_acids)
  three_di = read_fasta(family.three_di)
  if amino_acids.keys() != three_di.keys():
    only_aa = sorted(amino_acids.keys() - three_di.keys())
    only_3di = sorted(three_di.keys() - amino_acids.keys())
    raise AlignmentError(
      f"FASTA identifiers differ for {family.stem}; "
      f"AA-only={only_aa[:5]}, 3Di-only={only_3di[:5]}"
    )
  for omaid, sequence in amino_acids.items():
    structural = three_di[omaid]
    if "-" in sequence or "-" in structural:
      raise AlignmentError(f"Input sequences must be unaligned: {family.stem}/{omaid}")
    if len(sequence) != len(structural):
      raise AlignmentError(
        f"AA/3Di length mismatch for {family.stem}/{omaid}: "
        f"{len(sequence)} != {len(structural)}"
      )
  return amino_acids, three_di


def find_foldmason(command: str) -> str:
  if os.sep in command:
    candidate = Path(command)
    if candidate.is_file() and os.access(candidate, os.X_OK):
      return str(candidate)
  else:
    resolved = shutil.which(command)
    if resolved:
      return resolved
  raise AlignmentError(f"FoldMason executable not found: {command}")


def run_command(
  command: list[str],
  started_at: float,
  label: str,
) -> None:
  try:
    process = subprocess.Popen(command)
    while True:
      try:
        return_code = process.wait(timeout=60)
        break
      except subprocess.TimeoutExpired:
        progress(started_at, f"{label}: still running")
  except OSError as error:
    raise AlignmentError(f"Could not run {command[0]}: {error}") from error
  if return_code:
    raise CommandError(command, return_code)


def output_paths(output_dir: Path, stem: str) -> tuple[Path, Path, Path]:
  prefix = output_dir / "alignments" / stem
  return Path(f"{prefix}_aa.fa"), Path(f"{prefix}_3di.fa"), Path(f"{prefix}.nw")


def validate_alignment(
  family: FamilyInput,
  amino_acids: dict[str, str],
  three_di: dict[str, str],
  aa_path: Path,
  three_di_path: Path,
) -> None:
  aligned_aa = read_fasta(aa_path)
  aligned_3di = read_fasta(three_di_path)
  expected_ids = amino_acids.keys()
  if aligned_aa.keys() != expected_ids or aligned_3di.keys() != expected_ids:
    raise AlignmentError(f"FoldMason output identifiers differ from OMA for {family.stem}")
  lengths = {len(sequence) for sequence in aligned_aa.values()}
  if len(lengths) != 1:
    raise AlignmentError(f"FoldMason amino-acid output is not rectangular for {family.stem}")
  alignment_length = next(iter(lengths))
  if any(len(sequence) != alignment_length for sequence in aligned_3di.values()):
    raise AlignmentError(f"FoldMason AA and 3Di alignment lengths differ for {family.stem}")
  for omaid in amino_acids:
    if aligned_aa[omaid].replace("-", "").upper() != amino_acids[omaid].upper():
      raise AlignmentError(f"FoldMason changed the amino-acid sequence for {family.stem}/{omaid}")
    if aligned_3di[omaid].replace("-", "").upper() != three_di[omaid].upper():
      raise AlignmentError(f"FoldMason changed the 3Di sequence for {family.stem}/{omaid}")


def align_family(
  foldmason: str,
  family: FamilyInput,
  amino_acids: dict[str, str],
  three_di: dict[str, str],
  output_dir: Path,
  threads: int,
  refine_iters: int,
  keep_work: bool,
  started_at: float,
) -> tuple[Path, Path, Path]:
  work_root = output_dir / ".work"
  work_root.mkdir(parents=True, exist_ok=True)
  work_dir = Path(tempfile.mkdtemp(prefix=f"{family.stem}.", dir=work_root))
  database = work_dir / "family_db"
  staged_prefix = work_dir / "alignment"
  try:
    common = ["--shuffle", "0", "--threads", str(threads), "-v", "1"]
    run_command(
      [foldmason, "base:createdb", str(family.amino_acids), str(database), *common],
      started_at,
      f"{family.stem}: building amino-acid database",
    )
    run_command(
      [foldmason, "base:createdb", str(family.three_di), f"{database}_ss", *common],
      started_at,
      f"{family.stem}: building 3Di database",
    )
    command = [
      foldmason,
      "structuremsa",
      str(database),
      str(staged_prefix),
      "--threads",
      str(threads),
      "--refine-iters",
      str(refine_iters),
      "-v",
      "1",
    ]
    try:
      run_command(command, started_at, f"{family.stem}: FoldMason alignment")
    except CommandError as error:
      if error.return_code != -signal.SIGSEGV or threads == 1:
        raise
      progress(
        started_at,
        f"{family.stem}: FoldMason crashed with {threads} threads; retrying with 1",
      )
      retry_command = command.copy()
      retry_command[retry_command.index("--threads") + 1] = "1"
      run_command(retry_command, started_at, f"{family.stem}: FoldMason single-thread retry")
    staged = (
      Path(f"{staged_prefix}_aa.fa"),
      Path(f"{staged_prefix}_3di.fa"),
      Path(f"{staged_prefix}.nw"),
    )
    if any(not path.is_file() for path in staged):
      raise AlignmentError(f"FoldMason did not create all expected outputs for {family.stem}")
    validate_alignment(family, amino_acids, three_di, staged[0], staged[1])
    final = output_paths(output_dir, family.stem)
    for source, destination in zip(staged, final):
      atomic_write(destination, source.read_bytes())
  except BaseException:
    print(f"Retained failed work directory: {work_dir}", file=sys.stderr, flush=True)
    raise
  if not keep_work:
    shutil.rmtree(work_dir)
  return final


def tsv_text(
  rows: Iterable[dict[str, object]],
  fields: Sequence[str] = MANIFEST_FIELDS,
) -> str:
  from io import StringIO

  stream = StringIO(newline="")
  writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", extrasaction="ignore")
  writer.writeheader()
  writer.writerows(rows)
  return stream.getvalue()


def build_parser() -> argparse.ArgumentParser:
  parser = argparse.ArgumentParser(
    description="Align paired OMA amino-acid and ProstT5 3Di sequences with FoldMason."
  )
  parser.add_argument("--input-dir", required=True, type=Path)
  parser.add_argument(
    "--output-dir",
    type=Path,
    help="default: INPUT_DIR/structure-alignments",
  )
  parser.add_argument("--family", action="append", default=[], help="family stem; repeatable")
  parser.add_argument("--max-families", type=int)
  parser.add_argument("--threads", type=int, default=os.cpu_count() or 1)
  parser.add_argument("--foldmason", default="foldmason")
  parser.add_argument("--refine-iters", type=int, default=0)
  parser.add_argument("--force", action="store_true", help="replace completed alignments")
  parser.add_argument("--keep-work", action="store_true", help="keep intermediate Foldseek databases")
  parser.add_argument(
    "--stop-on-error",
    action="store_true",
    help="stop after the first failed family instead of continuing",
  )
  return parser


def validate_args(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
  if args.threads < 1:
    parser.error("--threads must be positive")
  if args.refine_iters < 0:
    parser.error("--refine-iters cannot be negative")
  if args.max_families is not None and args.max_families < 1:
    parser.error("--max-families must be positive")


def main(argv: Sequence[str] | None = None) -> int:
  started_at = time.monotonic()
  parser = build_parser()
  args = parser.parse_args(argv)
  validate_args(parser, args)
  input_dir = args.input_dir.resolve()
  output_dir = (
    args.output_dir.resolve()
    if args.output_dir is not None
    else input_dir / "structure-alignments"
  )
  try:
    families = discover_families(input_dir)
    if args.family:
      requested = {safe_name(value) for value in args.family}
      families = [family for family in families if family.stem in requested]
      missing = requested - {family.stem for family in families}
      if missing:
        raise AlignmentError(f"Requested families not found: {', '.join(sorted(missing))}")
    if args.max_families is not None:
      families = families[:args.max_families]
    if not families:
      raise AlignmentError("No families selected")
    foldmason = find_foldmason(args.foldmason)
  except AlignmentError as error:
    parser.error(str(error))

  output_dir.mkdir(parents=True, exist_ok=True)
  manifest_path = output_dir / "structure-alignments.tsv"
  errors_path = output_dir / "structure-alignment-errors.tsv"
  progress(started_at, f"Selected {len(families)} OMA gene families")
  manifest: list[dict[str, object]] = []
  errors: list[dict[str, object]] = []
  atomic_write_text(errors_path, tsv_text(errors, ERROR_FIELDS))
  failures = 0
  for index, family in enumerate(families, start=1):
    progress(started_at, f"[{index}/{len(families)}] {family.stem}")
    member_count: int | str = ""
    try:
      amino_acids, three_di = validate_family(family)
      member_count = len(amino_acids)
      aa_path, three_di_path, tree_path = output_paths(output_dir, family.stem)
      status = "aligned"
      if not args.force and aa_path.is_file() and three_di_path.is_file() and tree_path.is_file():
        validate_alignment(family, amino_acids, three_di, aa_path, three_di_path)
        status = "existing"
      else:
        progress(started_at, f"{family.stem}: aligning {len(amino_acids)} paired AA/3Di sequences")
        aa_path, three_di_path, tree_path = align_family(
          foldmason,
          family,
          amino_acids,
          three_di,
          output_dir,
          args.threads,
          args.refine_iters,
          args.keep_work,
          started_at,
        )
      manifest.append(
        {
          "family": family.stem,
          "member_count": len(amino_acids),
          "status": status,
          "amino_acid_alignment": str(aa_path.relative_to(output_dir)),
          "three_di_alignment": str(three_di_path.relative_to(output_dir)),
          "guide_tree": str(tree_path.relative_to(output_dir)),
        }
      )
      progress(started_at, f"{family.stem}: {status} {len(amino_acids)} sequences")
    except (AlignmentError, OSError, ValueError) as error:
      failures += 1
      manifest.append(
        {"family": family.stem, "member_count": member_count, "status": "error"}
      )
      errors.append(
        {"family": family.stem, "member_count": member_count, "error": str(error)}
      )
      atomic_write_text(errors_path, tsv_text(errors, ERROR_FIELDS))
      progress(started_at, f"{family.stem}: error; details: {errors_path}")
      if args.stop_on_error:
        atomic_write_text(manifest_path, tsv_text(manifest))
        break
    atomic_write_text(manifest_path, tsv_text(manifest))

  progress(
    started_at,
    f"Finished {len(manifest) - failures} families; {failures} failed; output: {output_dir}",
  )
  return 1 if failures else 0


if __name__ == "__main__":
  raise SystemExit(main())
