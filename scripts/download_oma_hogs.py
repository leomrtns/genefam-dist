#!/usr/bin/env python3
"""Download level-specific OMA HOGs as gene-family FASTA files.

The fundamental unit is an ancestral gene: a HOG induced at ``--level``.
Several such genes can belong to the same root HOG when a duplication predates
the selected ancestor.  ``--grouping root-hog`` merges those components while
retaining their level-specific HOG IDs in the member metadata.

Only the Python standard library is required.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import sys
import tempfile
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Iterator, Sequence
from urllib.error import HTTPError, URLError
from urllib.parse import quote, urlencode
from urllib.request import Request, urlopen


DEFAULT_API = "https://omabrowser.org/api"
TRANSIENT_HTTP_CODES = {429, 500, 502, 503, 504}
SEQUENCE_TYPES = ("protein", "cdna", "3di")


class OMAError(RuntimeError):
    """An OMA response was unavailable or did not match the expected schema."""


def _atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as handle:
            handle.write(text)
        os.replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


class OMAClient:
    """Small, cached and retrying client for the public OMA REST API."""

    def __init__(
        self,
        api_base: str = DEFAULT_API,
        cache_dir: Path | None = None,
        refresh_cache: bool = False,
        timeout: float = 60.0,
        retries: int = 5,
        started_at: float | None = None,
    ) -> None:
        self.api_base = api_base.rstrip("/")
        self.cache_dir = cache_dir
        self.refresh_cache = refresh_cache
        self.timeout = timeout
        self.retries = retries
        self.started_at = time.monotonic() if started_at is None else started_at
        self._cache_lock = threading.Lock()

    def _url(self, path: str, params: dict[str, Any] | None = None) -> str:
        url = f"{self.api_base}/{path.lstrip('/')}"
        if params:
            url = f"{url}?{urlencode(params)}"
        return url

    def _cache_path(self, url: str) -> Path | None:
        if self.cache_dir is None:
            return None
        digest = hashlib.sha256(url.encode("utf-8")).hexdigest()
        return self.cache_dir / digest[:2] / f"{digest}.json"

    def get_json(
        self, path: str, params: dict[str, Any] | None = None
    ) -> Any:
        url = self._url(path, params)
        cache_path = self._cache_path(url)
        if cache_path is not None and cache_path.exists() and not self.refresh_cache:
            try:
                return json.loads(cache_path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                pass

        request = Request(
            url,
            headers={
                "Accept": "application/json",
                "User-Agent": "genefam-dist-oma-downloader/1.0",
            },
        )
        for attempt in range(self.retries):
            try:
                with urlopen(request, timeout=self.timeout) as response:
                    payload = response.read().decode("utf-8")
                data = json.loads(payload)
                if cache_path is not None:
                    with self._cache_lock:
                        _atomic_write_text(cache_path, payload)
                return data
            except HTTPError as error:
                retry = error.code in TRANSIENT_HTTP_CODES
                detail = error.read().decode("utf-8", errors="replace")[:500]
                failure = OMAError(f"OMA returned HTTP {error.code} for {url}: {detail}")
            except (URLError, TimeoutError, json.JSONDecodeError) as error:
                retry = True
                failure = OMAError(f"Could not retrieve {url}: {error}")
            if not retry or attempt + 1 == self.retries:
                raise failure
            print(
                f"[{elapsed_minutes(self.started_at):.1f} min] "
                f"OMA request failed; retrying {url} "
                f"({attempt + 2}/{self.retries}): {failure}",
                file=sys.stderr,
                flush=True,
            )
            time.sleep(2**attempt)
        raise AssertionError("unreachable")


def iter_hogs_at_level(
    client: OMAClient, level: str, page_size: int = 1000
) -> Iterator[dict[str, Any]]:
    """Yield all HOGs representing genes in an ancestral genome."""
    started_at = getattr(client, "started_at", time.monotonic())
    page = 1
    while True:
        print(
            f"[{elapsed_minutes(started_at):.1f} min] "
            f"retrieving HOG metadata page {page}",
            file=sys.stderr,
            flush=True,
        )
        data = client.get_json(
            "/hog/", {"level": level, "page": page, "per_page": page_size}
        )
        if not isinstance(data, list):
            raise OMAError("OMA's HOG-list response is not a JSON list")
        for hog in data:
            if not isinstance(hog, dict) or "hog_id" not in hog:
                raise OMAError("OMA returned a malformed HOG record")
            yield hog
        if len(data) < page_size:
            break
        page += 1


def hogs_for_id(client: OMAClient, hog_id: str, level: str) -> list[dict[str, Any]]:
    """Resolve an ID to the HOG or induced HOGs present at ``level``."""
    encoded = quote(hog_id, safe=":.")
    data = client.get_json(f"/hog/{encoded}/", {"level": level})
    if isinstance(data, dict):
        data = [data]
    if not isinstance(data, list) or any(not isinstance(item, dict) for item in data):
        raise OMAError(f"Malformed HOG response for {hog_id}")
    return data


def get_members(client: OMAClient, hog_id: str, level: str) -> list[dict[str, Any]]:
    encoded = quote(hog_id, safe=":.")
    data = client.get_json(f"/hog/{encoded}/members/", {"level": level})
    if not isinstance(data, dict) or not isinstance(data.get("members"), list):
        raise OMAError(f"Malformed member response for {hog_id} at {level}")
    return data["members"]


def get_protein(client: OMAClient, member: dict[str, Any]) -> dict[str, Any]:
    entry_id = member.get("entry_nr") or member.get("omaid")
    if entry_id is None:
        raise OMAError("A HOG member has neither entry_nr nor omaid")
    data = client.get_json(f"/protein/{quote(str(entry_id), safe='')}/")
    if not isinstance(data, dict) or not data.get("omaid"):
        raise OMAError(f"Malformed protein response for {entry_id}")
    return data


def root_hog_id(hog: dict[str, Any]) -> str:
    value = hog.get("roothog_id")
    if value is not None:
        return f"HOG:F{int(value):07d}"
    return str(hog["hog_id"]).split(".", 1)[0]


def safe_family_name(family_id: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", family_id).strip("._")


def wrap_sequence(sequence: str, width: int = 80) -> str:
    return "\n".join(sequence[pos : pos + width] for pos in range(0, len(sequence), width))


def elapsed_minutes(started_at: float) -> float:
    return (time.monotonic() - started_at) / 60


def sequence_from_protein(protein: dict[str, Any], sequence_type: str) -> str:
    if sequence_type == "3di":
        structure = protein.get("structure")
        value = structure.get("sequence_3di", "") if isinstance(structure, dict) else ""
    else:
        key = "sequence" if sequence_type == "protein" else "cdna"
        value = protein.get(key, "")
    return str(value).replace("\n", "").replace(" ", "")


def family_fasta_paths(
    output_dir: Path, stem: str, sequence_type: str
) -> dict[str, Path]:
    extensions = {"protein": "faa", "cdna": "fna", "3di": "3di.fasta"}
    selected = SEQUENCE_TYPES if sequence_type == "all" else (sequence_type,)
    return {
        kind: output_dir / "families" / f"{stem}.{extensions[kind]}"
        for kind in selected
    }


def expected_member_count(family: Family) -> int | None:
    reported = [hog.get("nr_genes") for hog in family.components]
    if any(value is None for value in reported):
        return None
    return sum(int(float(value)) for value in reported)


def member_species(member: dict[str, Any]) -> dict[str, Any]:
    species = member.get("species")
    return species if isinstance(species, dict) else {}


@dataclass(frozen=True)
class Family:
    family_id: str
    root_hog: str
    components: tuple[dict[str, Any], ...]


def group_hogs(hogs: Iterable[dict[str, Any]], grouping: str) -> Iterator[Family]:
    if grouping == "ancestral-gene":
        for hog in hogs:
            yield Family(str(hog["hog_id"]), root_hog_id(hog), (hog,))
        return

    grouped: dict[str, list[dict[str, Any]]] = {}
    order: list[str] = []
    for hog in hogs:
        root = root_hog_id(hog)
        if root not in grouped:
            grouped[root] = []
            order.append(root)
        grouped[root].append(hog)
    for root in order:
        yield Family(root, root, tuple(grouped[root]))


def _tsv_text(fieldnames: Sequence[str], rows: Iterable[dict[str, Any]]) -> str:
    from io import StringIO

    stream = StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue()


MEMBER_FIELDS = (
    "omaid",
    "canonical_id",
    "entry_nr",
    "species_code",
    "species_name",
    "taxon_id",
    "ancestor_hog_id",
    "root_hog_id",
    "member_hog_id",
    "sequence_length",
    "sequence_md5",
)


def download_family(
    client: OMAClient,
    family: Family,
    level: str,
    output_dir: Path,
    sequence_type: str,
    workers: int,
    min_members: int,
    min_species: int,
    max_members: int | None,
    force: bool,
    started_at: float | None = None,
) -> dict[str, Any] | None:
    """Download one output family, returning its manifest row or None if filtered."""
    if started_at is None:
        started_at = time.monotonic()
    stem = safe_family_name(family.family_id)
    fasta_paths = family_fasta_paths(output_dir, stem, sequence_type)
    members_path = output_dir / "members" / f"{stem}.members.tsv"

    path_fields = {
        "protein_fasta": fasta_paths.get("protein"),
        "cdna_fasta": fasta_paths.get("cdna"),
        "three_di_fasta": fasta_paths.get("3di"),
    }
    relative_paths = {
        field: str(path.relative_to(output_dir)) if path is not None else ""
        for field, path in path_fields.items()
    }
    if all(path.exists() for path in fasta_paths.values()) and members_path.exists() and not force:
        return {
            "family_id": family.family_id,
            "root_hog_id": family.root_hog,
            "component_hogs": ",".join(str(h["hog_id"]) for h in family.components),
            "status": "existing",
            "fasta": next(path for path in relative_paths.values() if path),
            **relative_paths,
            "members_tsv": str(members_path.relative_to(output_dir)),
        }

    assignments: list[tuple[dict[str, Any], str]] = []
    seen_entries: set[Any] = set()
    for component_index, component in enumerate(family.components, start=1):
        component_id = str(component["hog_id"])
        print(
            f"[{elapsed_minutes(started_at):.1f} min] {family.family_id}: "
            f"retrieving member list for component "
            f"{component_index}/{len(family.components)} ({component_id})",
            file=sys.stderr,
            flush=True,
        )
        for member in get_members(client, component_id, level):
            key = member.get("entry_nr", member.get("omaid"))
            if key in seen_entries:
                raise OMAError(
                    f"Protein {key} occurs in more than one component of {family.family_id}"
                )
            seen_entries.add(key)
            assignments.append((member, component_id))

    print(
        f"[{elapsed_minutes(started_at):.1f} min] {family.family_id}: "
        f"found {len(assignments)} member genes",
        file=sys.stderr,
        flush=True,
    )

    species_codes = {member_species(member).get("code") for member, _ in assignments}
    species_codes.discard(None)
    if (
        len(assignments) < min_members
        or len(species_codes) < min_species
        or (max_members is not None and len(assignments) > max_members)
    ):
        return None

    print(
        f"[{elapsed_minutes(started_at):.1f} min] {family.family_id}: "
        f"downloading {len(assignments)} protein JSON records with {workers} workers",
        file=sys.stderr,
        flush=True,
    )

    proteins: list[dict[str, Any] | None] = [None] * len(assignments)
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(get_protein, client, pair[0]): index
            for index, pair in enumerate(assignments)
        }
        report_every = max(1, len(assignments) // 10)
        for downloaded, future in enumerate(as_completed(futures), start=1):
            proteins[futures[future]] = future.result()
            if downloaded == 1 or downloaded == len(assignments) or downloaded % report_every == 0:
                print(
                    f"[{elapsed_minutes(started_at):.1f} min] {family.family_id}: "
                    f"downloaded {downloaded}/{len(assignments)} protein JSON records",
                    file=sys.stderr,
                    flush=True,
                )

    fasta_parts: dict[str, list[str]] = {kind: [] for kind in fasta_paths}
    metadata: list[dict[str, Any]] = []
    for protein, (member, component_id) in zip(proteins, assignments):
        assert protein is not None
        omaid = str(protein["omaid"])
        sequences = {kind: sequence_from_protein(protein, kind) for kind in fasta_paths}
        for kind, sequence in sequences.items():
            if not sequence:
                raise OMAError(f"No {kind} sequence returned for {omaid}")
            fasta_parts[kind].append(f">{omaid}\n{wrap_sequence(sequence)}\n")
        species = member_species(member)
        metadata_sequence = sequences.get("protein") or next(iter(sequences.values()))
        metadata.append(
            {
                "omaid": omaid,
                "canonical_id": protein.get("canonicalid", member.get("canonicalid", "")),
                "entry_nr": protein.get("entry_nr", member.get("entry_nr", "")),
                "species_code": species.get("code", ""),
                "species_name": species.get("species", ""),
                "taxon_id": species.get("taxon_id", ""),
                "ancestor_hog_id": component_id,
                "root_hog_id": family.root_hog,
                "member_hog_id": member.get("oma_hog_id", ""),
                "sequence_length": len(metadata_sequence),
                "sequence_md5": protein.get("sequence_md5", member.get("sequence_md5", "")),
            }
        )

    for kind, fasta_path in fasta_paths.items():
        _atomic_write_text(fasta_path, "".join(fasta_parts[kind]))
    _atomic_write_text(members_path, _tsv_text(MEMBER_FIELDS, metadata))
    completeness = [float(h.get("completeness_score", 0.0)) for h in family.components]
    return {
        "family_id": family.family_id,
        "root_hog_id": family.root_hog,
        "component_hogs": ",".join(str(h["hog_id"]) for h in family.components),
        "component_count": len(family.components),
        "member_count": len(assignments),
        "species_count": len(species_codes),
        "min_component_completeness": min(completeness, default=""),
        "max_component_completeness": max(completeness, default=""),
        "description": family.components[0].get("description", ""),
        "status": "downloaded",
        "fasta": next(path for path in relative_paths.values() if path),
        **relative_paths,
        "members_tsv": str(members_path.relative_to(output_dir)),
    }


MANIFEST_FIELDS = (
    "family_id",
    "root_hog_id",
    "component_hogs",
    "component_count",
    "member_count",
    "species_count",
    "min_component_completeness",
    "max_component_completeness",
    "description",
    "status",
    "fasta",
    "protein_fasta",
    "cdna_fasta",
    "three_di_fasta",
    "members_tsv",
    "error",
)


def family_metadata_filter(
    family: Family,
    grouping: str,
    min_completeness: float,
    min_members: int,
    max_members: int | None,
) -> bool:
    """Cheap prefilter which never discards components from a merged family."""
    completeness = [float(hog.get("completeness_score", 0.0)) for hog in family.components]
    if grouping == "ancestral-gene":
        quality = completeness[0]
    else:
        # Select a root family if at least one ancestral copy is well supported,
        # then keep every copy belonging to that root family.
        quality = max(completeness, default=0.0)
    if quality < min_completeness:
        return False

    reported = [hog.get("nr_genes") for hog in family.components]
    if any(value is None for value in reported):
        return True
    member_count = sum(int(float(value)) for value in reported)
    return member_count >= min_members and (
        max_members is None or member_count <= max_members
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Download OMA HOG sequences for an ancestral genome."
    )
    parser.add_argument("--level", default="o__Enterobacterales", help="OMA taxonomic level")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--grouping",
        choices=("ancestral-gene", "root-hog"),
        default="ancestral-gene",
        help="keep inferred ancestral genes separate, or merge ancient duplicates by root HOG",
    )
    parser.add_argument(
        "--hog-id",
        action="append",
        default=[],
        help="download only this HOG (repeatable); older IDs may induce several HOGs at --level",
    )
    scope = parser.add_mutually_exclusive_group()
    scope.add_argument("--max-families", type=int, help="stop after this many output families")
    scope.add_argument("--all", action="store_true", help="explicitly allow all matching families")
    parser.add_argument(
        "--min-completeness",
        type=float,
        default=0.2,
        help=(
            "minimum HOG completeness (in root-hog mode, at least one component "
            "must pass and all of that root's components are retained)"
        ),
    )
    parser.add_argument("--min-members", type=int, default=4)
    parser.add_argument("--min-species", type=int, default=4)
    parser.add_argument("--max-members", type=int)
    parser.add_argument(
        "--sequence-type",
        choices=("all", *SEQUENCE_TYPES),
        default="all",
        help="sequence output to write; all reuses each protein JSON for protein, DNA, and 3Di",
    )
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--page-size", type=int, default=1000)
    parser.add_argument("--api-base", default=DEFAULT_API)
    parser.add_argument("--timeout", type=float, default=60.0)
    parser.add_argument("--retries", type=int, default=5)
    parser.add_argument("--refresh-cache", action="store_true")
    parser.add_argument("--force", action="store_true", help="replace completed family files")
    parser.add_argument("--dry-run", action="store_true", help="list matching HOG metadata only")
    return parser


def validate_args(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    if not args.hog_id and args.max_families is None and not args.all:
        parser.error("choose --max-families N for a trial run or --all for the full download")
    if args.max_families is not None and args.max_families < 1:
        parser.error("--max-families must be positive")
    if not 0 <= args.min_completeness <= 1:
        parser.error("--min-completeness must be between 0 and 1")
    for name in ("min_members", "min_species", "workers", "page_size", "retries"):
        if getattr(args, name) < 1:
            parser.error(f"--{name.replace('_', '-')} must be positive")


def selected_hogs(client: OMAClient, args: argparse.Namespace) -> Iterator[dict[str, Any]]:
    source: Iterable[dict[str, Any]]
    if args.hog_id:
        source = (
            hog
            for requested in args.hog_id
            for hog in hogs_for_id(client, requested, args.level)
        )
    else:
        source = iter_hogs_at_level(client, args.level, args.page_size)
    seen: set[str] = set()
    for hog in source:
        hog_id = str(hog["hog_id"])
        if hog_id in seen:
            continue
        seen.add(hog_id)
        yield hog


def main(argv: Sequence[str] | None = None) -> int:
    started_at = time.monotonic()
    parser = build_parser()
    args = parser.parse_args(argv)
    validate_args(parser, args)
    output_dir = args.output_dir.resolve()
    cache_dir = output_dir / ".cache"
    client = OMAClient(
        args.api_base,
        cache_dir,
        args.refresh_cache,
        args.timeout,
        args.retries,
        started_at,
    )

    hogs = selected_hogs(client, args)
    families = [
        family
        for family in group_hogs(hogs, args.grouping)
        if family_metadata_filter(
            family,
            args.grouping,
            args.min_completeness,
            args.min_members,
            args.max_members,
        )
    ]
    maximum = args.max_families if args.max_families is not None else len(families)
    print(
        f"[{elapsed_minutes(started_at):.1f} min] Found {len(families)} metadata-matching "
        f"gene families; expecting up to {min(maximum, len(families))} output families",
        file=sys.stderr,
        flush=True,
    )
    if args.dry_run:
        for index, family in enumerate(families, start=1):
            print(
                f"{family.family_id}\t{family.root_hog}\t"
                f"{','.join(str(h['hog_id']) for h in family.components)}"
            )
            if args.max_families is not None and index >= args.max_families:
                break
        return 0

    output_dir.mkdir(parents=True, exist_ok=True)
    manifest: list[dict[str, Any]] = []
    failures = 0
    completed = 0
    for candidate_index, family in enumerate(families, start=1):
        attempted = completed + failures
        if args.max_families is not None and attempted >= args.max_families:
            break
        expected_genes = expected_member_count(family)
        expectation = "unknown number of" if expected_genes is None else str(expected_genes)
        print(
            f"[{elapsed_minutes(started_at):.1f} min] "
            f"[{attempted + 1}/{min(maximum, len(families))}; "
            f"candidate {candidate_index}/{len(families)}] {family.family_id}: "
            f"expecting {expectation} member genes",
            file=sys.stderr,
            flush=True,
        )
        try:
            row = download_family(
                client,
                family,
                args.level,
                output_dir,
                args.sequence_type,
                args.workers,
                args.min_members,
                args.min_species,
                args.max_members,
                args.force,
                started_at,
            )
            if row is None:
                continue
            manifest.append(row)
            completed += 1
            _atomic_write_text(
                output_dir / "families.tsv", _tsv_text(MANIFEST_FIELDS, manifest)
            )
        except (OMAError, OSError, ValueError) as error:
            failures += 1
            manifest.append(
                {
                    "family_id": family.family_id,
                    "root_hog_id": family.root_hog,
                    "component_hogs": ",".join(
                        str(h["hog_id"]) for h in family.components
                    ),
                    "status": "error",
                    "error": str(error),
                }
            )
            print(f"error: {error}", file=sys.stderr)
            _atomic_write_text(
                output_dir / "families.tsv", _tsv_text(MANIFEST_FIELDS, manifest)
            )

    _atomic_write_text(output_dir / "families.tsv", _tsv_text(MANIFEST_FIELDS, manifest))
    run_info = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "api_base": args.api_base,
        "level": args.level,
        "grouping": args.grouping,
        "sequence_type": args.sequence_type,
        "min_completeness": args.min_completeness,
        "min_members": args.min_members,
        "min_species": args.min_species,
        "max_members": args.max_members,
        "requested_hog_ids": args.hog_id,
        "completed_families": completed,
        "failed_families": failures,
        "candidate_families": len(families),
    }
    _atomic_write_text(output_dir / "run.json", json.dumps(run_info, indent=2) + "\n")
    print(
        f"[{elapsed_minutes(started_at):.1f} min] Wrote {completed} families to {output_dir}; "
        f"kept reusable JSON cache in {cache_dir}",
        file=sys.stderr,
    )
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
