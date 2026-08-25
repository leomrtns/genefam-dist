#!/usr/bin/env python3
"""Download OMA members' AlphaFold models and align each family with FoldMason.

The input directory is an output directory from ``download_oma_hogs.py``. Each
family must have ``families/<name>.faa`` and ``families/<name>.3di.fasta``.
Member tables provide canonical IDs when present; OMA protein JSON files in
``.cache`` are used as a fallback.

Only the Python standard library is required. FoldMason is an external runtime
dependency unless ``--download-only`` is used.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Sequence
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import Request, urlopen


DEFAULT_AFDB_API = "https://alphafold.ebi.ac.uk/api/prediction"
TRANSIENT_HTTP_CODES = {429, 500, 502, 503, 504}


class AlignmentError(RuntimeError):
  """An input, download, or FoldMason result was invalid."""


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


def elapsed_minutes(started_at: float) -> float:
  return (time.monotonic() - started_at) / 60


def progress(started_at: float, message: str) -> None:
  print(
    f"[{elapsed_minutes(started_at):.1f} min] {message}",
    file=sys.stderr,
    flush=True,
  )


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


def read_member_ids(path: Path) -> dict[str, str]:
  try:
    with path.open(encoding="utf-8", newline="") as handle:
      rows = csv.DictReader(handle, delimiter="\t")
      if rows.fieldnames is None or "omaid" not in rows.fieldnames:
        raise AlignmentError(f"Missing omaid column in {path}")
      return {
        str(row.get("omaid", "")): str(row.get("canonical_id", ""))
        for row in rows
        if row.get("omaid")
      }
  except OSError as error:
    raise AlignmentError(f"Could not read {path}: {error}") from error


def cached_oma_proteins(cache_dir: Path) -> dict[str, dict[str, Any]]:
  """Index cached OMA protein records by OMA ID."""
  proteins: dict[str, dict[str, Any]] = {}
  if not cache_dir.is_dir():
    return proteins
  for path in cache_dir.rglob("*.json"):
    try:
      data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
      continue
    if isinstance(data, dict) and data.get("omaid") and data.get("sequence"):
      proteins[str(data["omaid"])] = data
  return proteins


@dataclass(frozen=True)
class FamilyInput:
  stem: str
  amino_acids: Path
  three_di: Path
  members: Path | None


def discover_families(input_dir: Path) -> list[FamilyInput]:
  family_dir = input_dir / "families"
  member_dir = input_dir / "members"
  families: list[FamilyInput] = []
  for amino_acids in sorted(family_dir.glob("*.faa")):
    stem = amino_acids.name.removesuffix(".faa")
    three_di = family_dir / f"{stem}.3di.fasta"
    members = member_dir / f"{stem}.members.tsv"
    if not three_di.is_file():
      raise AlignmentError(f"Missing companion input for {amino_acids}: {three_di}")
    families.append(
      FamilyInput(stem, amino_acids, three_di, members if members.is_file() else None)
    )
  if not families:
    raise AlignmentError(f"No family .faa files found under {family_dir}")
  return families


def validate_family_sequences(
  family: FamilyInput,
) -> tuple[dict[str, str], dict[str, str]]:
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
    if "-" in sequence or "-" in three_di[omaid]:
      raise AlignmentError(f"Input sequences must be unaligned: {family.stem}/{omaid}")
    if len(sequence) != len(three_di[omaid]):
      raise AlignmentError(
        f"AA/3Di length mismatch for {family.stem}/{omaid}: "
        f"{len(sequence)} != {len(three_di[omaid])}"
      )
  return amino_acids, three_di


class AFDBClient:
  """Cached and retrying client for AlphaFold DB metadata and model files."""

  def __init__(
    self,
    api_base: str,
    cache_dir: Path,
    timeout: float,
    retries: int,
    refresh: bool,
    started_at: float,
  ) -> None:
    self.api_base = api_base.rstrip("/")
    self.cache_dir = cache_dir
    self.timeout = timeout
    self.retries = retries
    self.refresh = refresh
    self.started_at = started_at
    self.cache_lock = threading.Lock()

  def request(self, url: str, accept: str) -> bytes:
    request = Request(
      url,
      headers={
        "Accept": accept,
        "User-Agent": "genefam-dist-structure-aligner/1.0",
      },
    )
    for attempt in range(self.retries):
      try:
        with urlopen(request, timeout=self.timeout) as response:
          return response.read()
      except HTTPError as error:
        retry = error.code in TRANSIENT_HTTP_CODES
        failure = AlignmentError(f"HTTP {error.code} for {url}")
      except (URLError, TimeoutError, OSError) as error:
        retry = True
        failure = AlignmentError(f"Could not retrieve {url}: {error}")
      if not retry or attempt + 1 == self.retries:
        raise failure
      progress(
        self.started_at,
        f"request failed; retrying {url} ({attempt + 2}/{self.retries}): {failure}",
      )
      time.sleep(2**attempt)
    raise AssertionError("unreachable")

  def metadata(self, accession: str) -> list[dict[str, Any]]:
    digest = hashlib.sha256(accession.encode("utf-8")).hexdigest()
    path = self.cache_dir / "metadata" / f"{digest}.json"
    if path.is_file() and not self.refresh:
      try:
        data = json.loads(path.read_text(encoding="utf-8"))
        if isinstance(data, list):
          return data
      except (OSError, json.JSONDecodeError):
        pass
    url = f"{self.api_base}/{quote(accession, safe='')}"
    payload = self.request(url, "application/json")
    try:
      data = json.loads(payload.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
      raise AlignmentError(f"Invalid AlphaFold metadata for {accession}: {error}") from error
    if not isinstance(data, list):
      raise AlignmentError(f"Unexpected AlphaFold metadata for {accession}")
    with self.cache_lock:
      atomic_write(path, payload)
    return data

  def model(self, url: str, path: Path) -> None:
    if path.is_file() and path.stat().st_size > 0 and not self.refresh:
      return
    payload = self.request(url, "chemical/x-mmcif,text/plain")
    if not payload.lstrip().startswith(b"data_"):
      raise AlignmentError(f"Downloaded AlphaFold model is not mmCIF: {url}")
    with self.cache_lock:
      atomic_write(path, payload)


def select_model(
  records: Iterable[dict[str, Any]],
  accession: str,
  expected_sequence: str,
) -> dict[str, Any]:
  exact = [
    record
    for record in records
    if record.get("cifUrl")
    and str(record.get("sequence") or record.get("uniprotSequence") or "")
    == expected_sequence
  ]
  if not exact:
    raise AlignmentError(
      f"AlphaFold DB has no model whose sequence exactly matches OMA {accession}"
    )
  exact.sort(
    key=lambda record: (
      bool(record.get("isComplex")),
      -int(record.get("latestVersion") or 0),
      str(record.get("modelEntityId") or ""),
    )
  )
  return exact[0]


def canonical_ids(
  family: FamilyInput,
  amino_acids: dict[str, str],
  member_identifiers: dict[str, str],
  cached_proteins: dict[str, dict[str, Any]],
) -> dict[str, str]:
  result: dict[str, str] = {}
  for omaid in amino_acids:
    canonical = member_identifiers.get(omaid, "")
    if not canonical:
      canonical = str(cached_proteins.get(omaid, {}).get("canonicalid", ""))
    if not canonical:
      raise AlignmentError(
        f"No canonical/UniProt ID for {family.stem}/{omaid} in member table or OMA cache"
      )
    result[omaid] = canonical
  return result


def download_family_models(
  client: AFDBClient,
  family: FamilyInput,
  amino_acids: dict[str, str],
  identifiers: dict[str, str],
  structure_dir: Path,
  workers: int,
  started_at: float,
) -> list[Path]:
  structure_dir.mkdir(parents=True, exist_ok=True)

  def download_one(omaid: str) -> Path:
    path = structure_dir / f"{safe_name(omaid)}.cif"
    accession = identifiers[omaid]
    record = select_model(client.metadata(accession), accession, amino_acids[omaid])
    client.model(str(record["cifUrl"]), path)
    return path

  paths: list[Path] = []
  with ThreadPoolExecutor(max_workers=workers) as executor:
    futures = {executor.submit(download_one, omaid): omaid for omaid in amino_acids}
    total = len(futures)
    report_every = max(1, total // 10)
    for completed, future in enumerate(as_completed(futures), start=1):
      omaid = futures[future]
      try:
        paths.append(future.result())
      except Exception as error:
        raise AlignmentError(f"Could not download model for {family.stem}/{omaid}: {error}") \
          from error
      if completed == 1 or completed == total or completed % report_every == 0:
        progress(
          started_at,
          f"{family.stem}: prepared {completed}/{total} AlphaFold models",
        )
  return sorted(paths)


def find_foldmason(command: str) -> str:
  if os.sep in command:
    path = Path(command)
    if path.is_file() and os.access(path, os.X_OK):
      return str(path)
  else:
    resolved = shutil.which(command)
    if resolved:
      return resolved
  raise AlignmentError(
    f"FoldMason executable not found: {command}; install it or use --download-only"
  )


def validate_foldmason_output(
  family: FamilyInput,
  output_prefix: Path,
  expected_amino_acids: dict[str, str],
) -> None:
  aa_path = Path(f"{output_prefix}_aa.fa")
  three_di_path = Path(f"{output_prefix}_3di.fa")
  if not aa_path.is_file() or not three_di_path.is_file():
    raise AlignmentError(
      f"FoldMason did not create expected outputs {aa_path.name} and {three_di_path.name}"
    )
  aligned_aa = read_fasta(aa_path)
  aligned_3di = read_fasta(three_di_path)
  if aligned_aa.keys() != aligned_3di.keys():
    raise AlignmentError(f"FoldMason AA and 3Di identifiers differ for {family.stem}")
  if aligned_aa.keys() != expected_amino_acids.keys():
    raise AlignmentError(f"FoldMason output identifiers differ from OMA for {family.stem}")
  lengths = {len(sequence) for sequence in aligned_aa.values()}
  if len(lengths) != 1 or any(len(aligned_3di[key]) not in lengths for key in aligned_3di):
    raise AlignmentError(f"FoldMason outputs are not rectangular for {family.stem}")
  for omaid, sequence in aligned_aa.items():
    if sequence.replace("-", "").upper() != expected_amino_acids[omaid].upper():
      raise AlignmentError(f"FoldMason changed the amino-acid sequence for {family.stem}/{omaid}")


def run_foldmason(
  foldmason: str,
  family: FamilyInput,
  structure_dir: Path,
  output_prefix: Path,
  expected_amino_acids: dict[str, str],
  report_mode: int,
) -> None:
  output_prefix.parent.mkdir(parents=True, exist_ok=True)
  temporary_root = output_prefix.parent.parent / ".foldmason"
  temporary_root.mkdir(parents=True, exist_ok=True)
  with tempfile.TemporaryDirectory(prefix=f"{family.stem}.", dir=temporary_root) as temporary:
    work_dir = Path(temporary)
    staged_prefix = work_dir / "result"
    foldmason_tmp = work_dir / "tmp"
    command = [
      foldmason,
      "easy-msa",
      str(structure_dir),
      str(staged_prefix),
      str(foldmason_tmp),
      "--report-mode",
      str(report_mode),
    ]
    try:
      subprocess.run(command, check=True)
    except subprocess.CalledProcessError as error:
      raise AlignmentError(
        f"FoldMason failed for {family.stem} with exit status {error.returncode}"
      ) from error
    validate_foldmason_output(family, staged_prefix, expected_amino_acids)
    staged_outputs = [path for path in work_dir.glob("result*") if path.is_file()]
    if not staged_outputs:
      raise AlignmentError(f"FoldMason created no output files for {family.stem}")
    for staged in staged_outputs:
      suffix = staged.name.removeprefix("result")
      os.replace(staged, Path(f"{output_prefix}{suffix}"))


MANIFEST_FIELDS = (
  "family",
  "member_count",
  "status",
  "structures",
  "amino_acid_alignment",
  "three_di_alignment",
  "error",
)


def tsv_text(rows: Iterable[dict[str, Any]]) -> str:
  from io import StringIO

  stream = StringIO(newline="")
  writer = csv.DictWriter(
    stream,
    fieldnames=MANIFEST_FIELDS,
    delimiter="\t",
    extrasaction="ignore",
  )
  writer.writeheader()
  writer.writerows(rows)
  return stream.getvalue()


def build_parser() -> argparse.ArgumentParser:
  parser = argparse.ArgumentParser(
    description="Download AlphaFold models for OMA families and align them with FoldMason."
  )
  parser.add_argument("--input-dir", required=True, type=Path)
  parser.add_argument(
    "--output-dir",
    type=Path,
    help="default: INPUT_DIR/structure-alignments",
  )
  parser.add_argument(
    "--family",
    action="append",
    default=[],
    help="family ID or filename stem (repeatable)",
  )
  parser.add_argument("--max-families", type=int)
  parser.add_argument("--workers", type=int, default=4)
  parser.add_argument("--foldmason", default="foldmason")
  parser.add_argument("--afdb-api", default=DEFAULT_AFDB_API)
  parser.add_argument("--timeout", type=float, default=60.0)
  parser.add_argument("--retries", type=int, default=5)
  parser.add_argument("--refresh", action="store_true", help="redownload models and realign")
  parser.add_argument(
    "--download-only",
    action="store_true",
    help="download and validate structures without running FoldMason",
  )
  parser.add_argument(
    "--report-mode",
    choices=(0, 1, 2),
    type=int,
    default=1,
    help="FoldMason report mode (default: 1, interactive HTML)",
  )
  return parser


def validate_args(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
  if args.workers < 1:
    parser.error("--workers must be positive")
  if args.retries < 1:
    parser.error("--retries must be positive")
  if args.timeout <= 0:
    parser.error("--timeout must be positive")
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
    foldmason = "" if args.download_only else find_foldmason(args.foldmason)
  except AlignmentError as error:
    parser.error(str(error))

  output_dir.mkdir(parents=True, exist_ok=True)
  cached_proteins: dict[str, dict[str, Any]] | None = None
  client = AFDBClient(
    args.afdb_api,
    output_dir / ".cache" / "afdb",
    args.timeout,
    args.retries,
    args.refresh,
    started_at,
  )
  progress(started_at, f"Selected {len(families)} OMA gene families")
  manifest: list[dict[str, Any]] = []
  failures = 0
  for index, family in enumerate(families, start=1):
    progress(started_at, f"[{index}/{len(families)}] {family.stem}")
    structure_dir = output_dir / "structures" / family.stem
    output_prefix = output_dir / "alignments" / family.stem
    aa_alignment = Path(f"{output_prefix}_aa.fa")
    three_di_alignment = Path(f"{output_prefix}_3di.fa")
    try:
      amino_acids, _three_di = validate_family_sequences(family)
      member_identifiers = (
        read_member_ids(family.members) if family.members is not None else {}
      )
      if any(not member_identifiers.get(omaid) for omaid in amino_acids):
        if cached_proteins is None:
          progress(started_at, "indexing cached OMA protein JSON for missing canonical IDs")
          cached_proteins = cached_oma_proteins(input_dir / ".cache")
      identifiers = canonical_ids(
        family,
        amino_acids,
        member_identifiers,
        cached_proteins or {},
      )
      model_paths = download_family_models(
        client,
        family,
        amino_acids,
        identifiers,
        structure_dir,
        args.workers,
        started_at,
      )
      status = "downloaded"
      if not args.download_only:
        if aa_alignment.is_file() and three_di_alignment.is_file() and not args.refresh:
          validate_foldmason_output(family, output_prefix, amino_acids)
          status = "existing"
        else:
          progress(started_at, f"{family.stem}: running FoldMason")
          run_foldmason(
            foldmason,
            family,
            structure_dir,
            output_prefix,
            amino_acids,
            args.report_mode,
          )
          status = "aligned"
      manifest.append(
        {
          "family": family.stem,
          "member_count": len(amino_acids),
          "status": status,
          "structures": str(structure_dir.relative_to(output_dir)),
          "amino_acid_alignment": (
            str(aa_alignment.relative_to(output_dir)) if aa_alignment.is_file() else ""
          ),
          "three_di_alignment": (
            str(three_di_alignment.relative_to(output_dir))
            if three_di_alignment.is_file()
            else ""
          ),
        }
      )
      progress(
        started_at,
        f"{family.stem}: {status} ({len(model_paths)} models)",
      )
    except (AlignmentError, OSError, ValueError) as error:
      failures += 1
      manifest.append(
        {
          "family": family.stem,
          "status": "error",
          "structures": str(structure_dir.relative_to(output_dir)),
          "error": str(error),
        }
      )
      progress(started_at, f"{family.stem}: error: {error}")
    atomic_write_text(output_dir / "structure-alignments.tsv", tsv_text(manifest))

  progress(
    started_at,
    f"Finished {len(families) - failures} families; {failures} failed; output: {output_dir}",
  )
  return 1 if failures else 0


if __name__ == "__main__":
  raise SystemExit(main())
