"""One workflow for comparing gene-family files with reusable reference trees.

The native objects own their tree pointers; Newick is only used at input/output
boundaries. This module needs neither NumPy nor DendroPy.
"""

import csv
import glob
import math
import os
from pathlib import Path
import re
import tempfile

from . import _treesignalc


LENGTH_METRICS = ("patristic", "wrf", "kf")
TOPOLOGY_METRICS = ("duplication", "loss", "deep_coalescence", "rf", "hdist", "spr")


def _tree_space(trees, rooted=False):
    if isinstance(trees, str):
        return _treesignalc.TreeSpace(trees, rooted=rooted)
    if not isinstance(trees, _treesignalc.TreeSpace):
        raise TypeError("trees must be Newick text or a TreeSpace")
    if not trees.ntrees:
        raise ValueError("TreeSpace is uninitialized")
    return trees


def _metrics(metrics):
    if isinstance(metrics, str):
        metrics = (metrics,)
    metrics = tuple(metrics)
    if not metrics:
        raise ValueError("select at least one metric")
    if len(set(metrics)) != len(metrics):
        raise ValueError("metrics must not contain duplicates")
    unknown = set(metrics).difference((*LENGTH_METRICS, "topology"))
    if unknown:
        raise ValueError("unknown metrics: " + ", ".join(sorted(unknown)))
    return metrics


def _rms_difference(first, second, normalize):
    if normalize:
        first_mean = math.fsum(value / len(first) for value in first)
        second_mean = math.fsum(value / len(second) for value in second)
        first = [value / first_mean for value in first] if first_mean else first
        second = [value / second_mean for value in second] if second_mean else second
    scale = math.sqrt(len(first))
    return math.hypot(*((x - y) / scale for x, y in zip(first, second)))


class TreeAnalysis:
    """Retain reference trees once and compare individual gene-family trees.

    ``compare`` returns one tuple per selected metric, in reference-tree order.
    Patristic distances use the species shared by the family and references.
    With paralogs, supply an explicit ``gene_to_species`` mapping and choose
    minimum (default) or mean cross-copy distances. Weighted RF and KF compare
    single-copy, identically labelled taxa on the induced common leaf set.
    """

    def __init__(self, references):
        self.references = _tree_space(references, rooted=True)

    def compare(self, gene, metrics=LENGTH_METRICS, normalize=False,
                gene_to_species=None, aggregation="min"):
        metrics = _metrics(metrics)
        if aggregation not in {"min", "mean"}:
            raise ValueError("aggregation must be 'min' or 'mean'")
        gene = _tree_space(gene)
        if gene.ntrees != 1:
            raise ValueError("each gene family must contain exactly one tree")
        if gene_to_species is not None:
            missing = set(gene.taxa).difference(gene_to_species)
            if missing:
                raise ValueError("missing gene-to-species mapping for: " + ", ".join(sorted(missing)))
            if any(not isinstance(gene_to_species[label], str) or not gene_to_species[label]
                   for label in gene.taxa):
                raise ValueError("mapped species labels must be nonempty strings")
            if any(metric != "patristic" for metric in metrics) and any(
                    gene_to_species[label] != label for label in gene.taxa):
                raise ValueError("gene-to-species mapping is supported for patristic distances; "
                                 "use metrics=('patristic',) for paralogous or renamed leaves")
        result = {}
        if "patristic" in metrics:
            result["patristic"] = self._patristic(gene, normalize, gene_to_species, aggregation)
        if "wrf" in metrics or "kf" in metrics:
            values = _treesignalc.from_tree_spaces_lengths(gene, self.references, normalize=bool(normalize))
            for metric, offset in (("wrf", 1), ("kf", 2)):
                if metric in metrics:
                    result[metric] = tuple(values[offset::3])
        if "topology" in metrics:
            function = (_treesignalc.from_tree_spaces_rescale if normalize
                        else _treesignalc.from_tree_spaces)
            result["topology"] = function(gene, self.references)
        return {metric: result[metric] for metric in metrics}

    def _patristic(self, gene, normalize, gene_to_species, aggregation):
        if not gene.has_branch_lengths or not self.references.has_branch_lengths:
            raise ValueError("patristic distances require complete nonnegative branch lengths "
                             "on the gene and reference trees")
        groups = {}
        for index, label in enumerate(gene.taxa):
            species = label if gene_to_species is None else gene_to_species[label]
            groups.setdefault(species, []).append(index)
        reference_ids = {label: index for index, label in enumerate(self.references.taxa)}
        common = sorted(set(groups).intersection(reference_ids))
        if len(common) < 2:
            raise ValueError("patristic distances require at least two shared species")
        pairs = [(a, b) for index, a in enumerate(common) for b in common[index + 1:]]
        gene_matrix = gene.patristic_matrix()
        projected = []
        for a, b in pairs:
            copies = [gene_matrix[i][j] for i in groups[a] for j in groups[b]]
            projected.append(min(copies) if aggregation == "min"
                             else math.fsum(value / len(copies) for value in copies))
        distances = []
        for index in range(self.references.ndistinct):
            matrix = self.references.patristic_matrix(index)
            reference = [matrix[reference_ids[a]][reference_ids[b]] for a, b in pairs]
            if not all(math.isfinite(value) for value in (*projected, *reference)):
                raise ValueError("patristic path distances exceed the finite numeric range")
            distances.append(_rms_difference(projected, reference, normalize))
        return tuple(distances)


def generate_references(species, count=20, n_spr=1, seed=0, branch_length=None):
    """Generate an SPR chain in the native library, with the supplied labels.

    ``count`` includes the initial random tree. A nonzero seed is reproducible.
    Edges have unit lengths (the artificial root splits one unit edge in half).
    ``branch_length`` rescales these synthetic edges; these are reference design
    choices, not estimates of species divergence times.
    """
    if isinstance(species, (str, bytes)):
        raise TypeError("species must be a sequence of labels, not a string")
    species = tuple(species)
    if len(species) < 4:
        raise ValueError("reference generation requires at least four species")
    if len(set(species)) != len(species):
        raise ValueError("species labels must be unique")
    if any(not isinstance(label, str) or re.fullmatch(r"[A-Za-z0-9_.|+\-/]+", label) is None
           for label in species):
        raise ValueError("species labels must be simple unquoted Newick labels")
    for value, name in ((count, "count"), (n_spr, "n_spr"), (seed, "seed")):
        if isinstance(value, bool) or not isinstance(value, int):
            raise TypeError(name + " must be an integer")
    if count < 1 or n_spr < 1 or not 0 <= seed <= 2**64 - 1:
        raise ValueError("count and n_spr must be positive; seed must be an unsigned 64-bit integer")
    if branch_length is not None:
        branch_length = float(branch_length)
        if not math.isfinite(branch_length) or branch_length < 0:
            raise ValueError("branch_length must be finite and nonnegative")
    generated = _treesignalc.generate_tree_space(len(species), count, n_spr, seed)
    mapping = {"T" + str(index + 1): label for index, label in enumerate(species)}
    newick = re.sub(r"(?<=[(,])T[0-9]+(?=[:),;])", lambda match: mapping[match.group()], generated.newick)
    if branch_length is not None:
        newick = re.sub(r":([^,);\s]+)", lambda match: ":" + format(float(match.group(1)) * branch_length, ".17g"), newick)
    return _treesignalc.TreeSpace(newick, rooted=True)


def _paths(paths):
    if isinstance(paths, (str, os.PathLike)):
        paths = [paths]
    expanded = []
    for pattern in paths:
        pattern = os.fspath(pattern)
        matches = sorted(glob.glob(pattern))
        if not matches:
            raise ValueError("no gene-tree files match: " + pattern)
        expanded.extend(Path(path) for path in matches)
    expanded = sorted(set(expanded), key=str)
    if not expanded:
        raise ValueError("at least one gene-tree file is required")
    names = [path.name for path in expanded]
    if len(set(names)) != len(names):
        raise ValueError("gene-tree filenames must be unique across input directories")
    return expanded


def _write_rows(output, rows):
    output = Path(output)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", newline="",
                                         dir=output.parent, prefix="." + output.name + ".", delete=False) as handle:
            temporary = Path(handle.name)
            writer = csv.DictWriter(handle, fieldnames=("family", "reference", "metric", "value"), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary, output)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def analyze_files(paths, references, output=None, metrics=LENGTH_METRICS,
                  normalize=False, gene_to_species=None, aggregation="min"):
    """Compare one Newick tree per file and optionally save a tidy TSV table.

    Each row has ``family`` (filename), 1-based ``reference``, ``metric``, and
    ``value`` fields. Inputs are sorted for repeatable output. The output is
    replaced only after all input trees have been successfully analyzed.
    """
    paths = _paths(paths)
    if output is not None and Path(output).resolve() in {path.resolve() for path in paths}:
        raise ValueError("output must not overwrite an input gene-tree file")
    analysis = references if isinstance(references, TreeAnalysis) else TreeAnalysis(references)
    rows = []
    for path in paths:
        try:
            distances = analysis.compare(path.read_text(encoding="utf-8"), metrics=metrics,
                                         normalize=normalize, gene_to_species=gene_to_species,
                                         aggregation=aggregation)
        except (OSError, ValueError, TypeError, _treesignalc.error) as error:
            raise ValueError(str(path) + ": " + str(error)) from error
        for metric, values in distances.items():
            for index, value in enumerate(values):
                if metric == "topology":
                    name = TOPOLOGY_METRICS[index % len(TOPOLOGY_METRICS)]
                    reference = index // len(TOPOLOGY_METRICS) + 1
                else:
                    name = metric
                    reference = index + 1
                rows.append({"family": path.name, "reference": reference,
                             "metric": name, "value": value})
    if output is not None:
        _write_rows(output, rows)
    return rows
