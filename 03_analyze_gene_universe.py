from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import sqlite3
import statistics
import sys
import tempfile
import urllib.error
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable

import numpy as np
from scipy.stats import mannwhitneyu


ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "Data"
LABEL_DIR = DATA_DIR / "processed_labels"
UNIVERSE_DIR = DATA_DIR / "gene_universe_analysis"
OUTPUT_DIR = DATA_DIR / "sl_analysis"
UNIVERSE_FILE = UNIVERSE_DIR / "candidate_gene_universe.tsv"
PAIR_SUMMARY_FILE = LABEL_DIR / "pair_summary_all.tsv"
BASELINE_FILE = LABEL_DIR / "pair_baseline_no_conflicts.tsv"
RAW_SL_FILE = DATA_DIR / "Human.SL.detailed.tsv"
RAW_NONSL_FILE = DATA_DIR / "Human.non.SL.detailed.tsv"
STRING_LINKS_FILE = DATA_DIR / "9606.protein.links.v12.5.txt.gz"
STRING_ALIASES_FILE = DATA_DIR / "9606.protein.aliases.v12.5.txt.gz"

GPROFILER_URL = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
GPROFILER_SOURCES = ["GO:BP", "GO:MF", "GO:CC", "REAC", "KEGG"]
STRING_GENE_ID_SOURCE = "Ensembl_HGNC_entrez_id"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Summarize SL labels and partners, run GO/pathway enrichment, "
            "and test pair-label association with STRING degree."
        )
    )
    parser.add_argument(
        "--string-min-score",
        type=int,
        default=700,
        help=(
            "Minimum STRING combined_score (0-1000) used to build the "
            "independent gene network. Default: 700."
        ),
    )
    parser.add_argument(
        "--skip-enrichment",
        action="store_true",
        help="Skip the online g:Profiler enrichment request.",
    )
    parser.add_argument(
        "--permutations",
        type=int,
        default=9999,
        help="Number of node-label permutations for the centrality test.",
    )
    return parser.parse_args()


def require_files(paths: Iterable[Path]) -> None:
    missing = [str(path) for path in paths if not path.is_file()]
    if missing:
        raise FileNotFoundError("Required input file(s) missing: " + ", ".join(missing))


def read_tsv(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as file:
        yield from csv.DictReader(file, delimiter="\t")


def write_tsv(path: Path, columns: list[str], rows: Iterable[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def canonical_pair(gene_a: str, gene_b: str) -> tuple[str, str] | None:
    gene_a = gene_a.strip()
    gene_b = gene_b.strip()
    if not gene_a or not gene_b or gene_a == gene_b:
        return None
    return tuple(sorted((gene_a, gene_b)))


def read_gene_universe(path: Path) -> tuple[dict[str, str], set[str]]:
    symbols: dict[str, str] = {}
    with path.open("r", encoding="utf-8-sig", newline="") as file:
        reader = csv.DictReader(file, delimiter="\t")
        required = {"entrez_id", "expression_column", "crispr_column"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path.name} missing columns: {sorted(missing)}")
        for row in reader:
            gene_id = row["entrez_id"].strip()
            if not gene_id:
                raise ValueError(f"Blank Entrez ID in {path}")
            symbols[gene_id] = (
                row["expression_column"].strip()
                or row["crispr_column"].strip()
                or gene_id
            )
    return symbols, set(symbols)


def read_raw_pair_file(path: Path) -> tuple[int, set[tuple[str, str]]]:
    row_count = 0
    pairs: set[tuple[str, str]] = set()
    for row in read_tsv(path):
        row_count += 1
        pair = canonical_pair(row.get("x:START_ID", ""), row.get("y:END_ID", ""))
        if pair:
            pairs.add(pair)
    return row_count, pairs


def read_labeled_pairs(path: Path) -> tuple[dict[tuple[str, str], int], Counter]:
    pairs: dict[tuple[str, str], int] = {}
    counts: Counter = Counter()
    for row in read_tsv(path):
        pair = canonical_pair(row.get("gene_a_id", ""), row.get("gene_b_id", ""))
        if not pair:
            continue
        status = row.get("label_status", "")
        counts[status] += 1
        if status in {"SL_ONLY", "NONSL_ONLY"}:
            pairs[pair] = 1 if status == "SL_ONLY" else 0
    return pairs, counts


def read_baseline_pairs(path: Path) -> dict[tuple[str, str], int]:
    pairs: dict[tuple[str, str], int] = {}
    for row in read_tsv(path):
        pair = canonical_pair(row.get("gene_a_id", ""), row.get("gene_b_id", ""))
        if not pair:
            continue
        label = row.get("label", "")
        if label not in {"0", "1"}:
            raise ValueError(f"Invalid baseline label {label!r} in {path}")
        pairs[pair] = int(label)
    return pairs


def partner_counts(
    pairs: Iterable[tuple[str, str]],
) -> dict[str, set[str]]:
    partners: dict[str, set[str]] = defaultdict(set)
    for gene_a, gene_b in pairs:
        partners[gene_a].add(gene_b)
        partners[gene_b].add(gene_a)
    return partners


def quantile_nearest_rank(sorted_values: list[int], probability: float) -> int:
    if not sorted_values:
        return 0
    index = max(0, math.ceil(probability * len(sorted_values)) - 1)
    return sorted_values[index]


def degree_summary(
    scope: str, degrees: dict[str, int], universe_size: int | None = None
) -> dict:
    values = list(degrees.values())
    if universe_size is not None:
        values.extend([0] * max(0, universe_size - len(values)))
    values.sort()
    nonzero = sum(value > 0 for value in values)
    return {
        "scope": scope,
        "gene_count_in_summary": len(values),
        "genes_with_at_least_one_partner": nonzero,
        "genes_with_zero_partners": len(values) - nonzero,
        "mean_partners": statistics.mean(values) if values else 0,
        "median_partners": statistics.median(values) if values else 0,
        "p90_nearest_rank": quantile_nearest_rank(values, 0.90),
        "p95_nearest_rank": quantile_nearest_rank(values, 0.95),
        "maximum_partners": max(values, default=0),
    }


def distribution_bins(
    scope: str, degrees: dict[str, int], universe_size: int | None = None
) -> list[dict]:
    values = list(degrees.values())
    if universe_size is not None:
        values.extend([0] * max(0, universe_size - len(values)))
    bins = [
        ("0", lambda value: value == 0),
        ("1", lambda value: value == 1),
        ("2-5", lambda value: 2 <= value <= 5),
        ("6-10", lambda value: 6 <= value <= 10),
        ("11-50", lambda value: 11 <= value <= 50),
        ("51-100", lambda value: 51 <= value <= 100),
        ("101-500", lambda value: 101 <= value <= 500),
        (">500", lambda value: value > 500),
    ]
    return [
        {
            "scope": scope,
            "partner_count_bin": name,
            "gene_count": sum(predicate(value) for value in values),
        }
        for name, predicate in bins
    ]


def write_distribution_plot(
    label_counts: Counter,
    degree_distributions: dict[str, dict[str, int]],
    candidate_count: int,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels = ["SL_ONLY", "NONSL_ONLY", "BOTH_LABELS"]
    values = [label_counts[label] for label in labels]
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    axes[0].bar(labels, values, color=["#c44e52", "#4c72b0", "#999999"])
    axes[0].set_title("Unique pair labels (all genes)")
    axes[0].set_ylabel("Number of unique pairs")
    axes[0].tick_params(axis="x", rotation=20)
    for index, value in enumerate(values):
        axes[0].text(index, value, f"{value:,}", ha="center", va="bottom")

    scopes = [
        ("Raw SL (observed genes)", degree_distributions["raw_sl"]),
        ("Baseline SL_ONLY (candidate universe)", degree_distributions["candidate_sl"]),
    ]
    for name, counts in scopes:
        positive_degrees = [
            degree for degree in counts.values() if degree > 0
        ]
        if positive_degrees:
            axes[1].hist(
                positive_degrees,
                bins=50,
                histtype="step",
                linewidth=1.8,
                label=name,
            )
    axes[1].set_xscale("log")
    axes[1].set_yscale("log")
    axes[1].set_title(
        f"SL partner counts (candidate genes: {candidate_count:,})"
    )
    axes[1].set_xlabel("Unique SL partners per gene (log scale)")
    axes[1].set_ylabel("Number of genes (log scale)")
    axes[1].legend()
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "sl_distribution.png", dpi=180)
    plt.close(fig)


def write_enrichment_query(
    query_genes: set[str], background_genes: set[str]
) -> None:
    write_tsv(
        OUTPUT_DIR / "go_query_gene_ids.tsv",
        ["entrez_id"],
        ({"entrez_id": gene_id} for gene_id in sorted(query_genes, key=int)),
    )
    write_tsv(
        OUTPUT_DIR / "go_background_gene_ids.tsv",
        ["entrez_id"],
        ({"entrez_id": gene_id} for gene_id in sorted(background_genes, key=int)),
    )


def run_gprofiler(
    query_genes: set[str], background_genes: set[str]
) -> tuple[int, str]:
    payload = {
        "organism": "hsapiens",
        "query": sorted(query_genes, key=int),
        "sources": GPROFILER_SOURCES,
        "user_threshold": 0.05,
        "significance_threshold_method": "fdr",
        "domain_scope": "custom",
        "background": sorted(background_genes, key=int),
        "no_evidences": True,
    }
    request = urllib.request.Request(
        GPROFILER_URL,
        data=json.dumps(payload).encode("utf-8"),
        headers={
            "Content-Type": "application/json",
            "User-Agent": "HUST-Project-I-SL-analysis/1.0",
        },
        method="POST",
    )
    try:
        with urllib.request.urlopen(request, timeout=120) as response:
            result = json.loads(response.read().decode("utf-8"))
    except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError) as error:
        error_text = f"{type(error).__name__}: {error}"
        (OUTPUT_DIR / "gprofiler_error.txt").write_text(
            error_text + "\n", encoding="utf-8"
        )
        return 0, error_text

    if "result" not in result:
        error_text = f"g:Profiler response has no 'result': {result}"
        (OUTPUT_DIR / "gprofiler_error.txt").write_text(
            error_text + "\n", encoding="utf-8"
        )
        return 0, error_text

    results = result["result"]
    columns = [
        "source",
        "native",
        "name",
        "p_value_fdr",
        "intersection_size",
        "query_size",
        "term_size",
        "effective_domain_size",
        "description",
    ]
    rows = []
    for item in results:
        rows.append(
            {
                "source": item.get("source", ""),
                "native": item.get("native", ""),
                "name": item.get("name", ""),
                "p_value_fdr": item.get("p_value", ""),
                "intersection_size": item.get("intersection_size", ""),
                "query_size": item.get("query_size", ""),
                "term_size": item.get("term_size", ""),
                "effective_domain_size": item.get("effective_domain_size", ""),
                "description": item.get("description", ""),
            }
        )
    write_tsv(OUTPUT_DIR / "go_pathway_enrichment.tsv", columns, rows)
    (OUTPUT_DIR / "gprofiler_request_metadata.json").write_text(
        json.dumps(
            {
                "endpoint": GPROFILER_URL,
                "organism": "hsapiens",
                "sources": GPROFILER_SOURCES,
                "query_size_submitted": len(query_genes),
                "background_size_submitted": len(background_genes),
                "domain_scope": "custom",
                "multiple_testing_method": "fdr",
                "significance_threshold": 0.05,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    error_path = OUTPUT_DIR / "gprofiler_error.txt"
    if error_path.exists():
        error_path.unlink()
    return len(rows), ""


def load_string_protein_mapping(candidate_genes: set[str]) -> tuple[dict, int]:
    protein_gene_candidates: dict[str, set[str]] = defaultdict(set)
    with gzip.open(STRING_ALIASES_FILE, "rt", encoding="utf-8") as file:
        reader = csv.DictReader(file, delimiter="\t")
        for row in reader:
            if row.get("source") != STRING_GENE_ID_SOURCE:
                continue
            gene_id = (row.get("alias") or "").strip()
            protein_id = (row.get("#string_protein_id") or "").strip()
            if gene_id in candidate_genes and protein_id:
                protein_gene_candidates[protein_id].add(gene_id)

    mapping = {
        protein_id: next(iter(gene_ids))
        for protein_id, gene_ids in protein_gene_candidates.items()
        if len(gene_ids) == 1
    }
    ambiguous = sum(
        len(gene_ids) > 1 for gene_ids in protein_gene_candidates.values()
    )
    return mapping, ambiguous


def build_string_gene_degrees(
    candidate_genes: set[str], minimum_score: int
) -> tuple[dict[str, int], dict[str, int]]:
    protein_to_gene, ambiguous_proteins = load_string_protein_mapping(
        candidate_genes
    )
    mapped_proteins = len(protein_to_gene)
    edge_count = 0

    with tempfile.TemporaryDirectory(
        prefix="sl_string_edges_", dir=OUTPUT_DIR
    ) as temporary_directory:
        database_path = Path(temporary_directory) / "edges.sqlite"
        connection = sqlite3.connect(database_path)
        try:
            connection.execute(
                "CREATE TABLE edge (gene_a TEXT NOT NULL, gene_b TEXT NOT NULL, "
                "PRIMARY KEY (gene_a, gene_b)) WITHOUT ROWID"
            )
            insert_sql = "INSERT OR IGNORE INTO edge (gene_a, gene_b) VALUES (?, ?)"
            batch = []

            with gzip.open(STRING_LINKS_FILE, "rt", encoding="utf-8") as file:
                next(file, None)
                for line in file:
                    fields = line.split()
                    if len(fields) != 3:
                        continue
                    protein_a, protein_b, score_text = fields
                    score = int(score_text)
                    if score < minimum_score:
                        continue
                    gene_a = protein_to_gene.get(protein_a)
                    gene_b = protein_to_gene.get(protein_b)
                    if not gene_a or not gene_b or gene_a == gene_b:
                        continue
                    batch.append(tuple(sorted((gene_a, gene_b))))
                    if len(batch) >= 50000:
                        connection.executemany(insert_sql, batch)
                        connection.commit()
                        batch.clear()

            if batch:
                connection.executemany(insert_sql, batch)
                connection.commit()

            edge_count = connection.execute(
                "SELECT COUNT(*) FROM edge"
            ).fetchone()[0]
            degree_counts: Counter = Counter()
            for gene_a, gene_b in connection.execute(
                "SELECT gene_a, gene_b FROM edge"
            ):
                degree_counts[gene_a] += 1
                degree_counts[gene_b] += 1
        finally:
            connection.close()

    degrees = {gene_id: degree_counts.get(gene_id, 0) for gene_id in candidate_genes}
    metadata = {
        "candidate_genes": len(candidate_genes),
        "mapped_string_proteins": mapped_proteins,
        "ambiguous_string_proteins_excluded": ambiguous_proteins,
        "gene_level_edges_after_threshold": edge_count,
        "string_minimum_combined_score": minimum_score,
        "alias_source": STRING_GENE_ID_SOURCE,
    }
    return degrees, metadata


def centrality_association(
    pairs: dict[tuple[str, str], int],
    candidate_genes: set[str],
    degree_by_gene: dict[str, int],
    permutation_count: int,
    random_seed: int = 20261109,
) -> tuple[list[dict], dict]:
    selected_pairs = [
        (pair, label)
        for pair, label in pairs.items()
        if pair[0] in candidate_genes and pair[1] in candidate_genes
    ]
    degree_by_label: dict[int, list[float]] = {0: [], 1: []}
    pair_rows = []
    node_count = len(candidate_genes)

    for (gene_a, gene_b), label in selected_pairs:
        degree_a = degree_by_gene.get(gene_a, 0)
        degree_b = degree_by_gene.get(gene_b, 0)
        mean_degree = (degree_a + degree_b) / 2
        mean_log_degree = (math.log1p(degree_a) + math.log1p(degree_b)) / 2
        degree_by_label[label].append(mean_degree)
        pair_rows.append(
            {
                "gene_a_id": gene_a,
                "gene_b_id": gene_b,
                "label": label,
                "gene_a_string_degree": degree_a,
                "gene_b_string_degree": degree_b,
                "mean_endpoint_degree": mean_degree,
                "max_endpoint_degree": max(degree_a, degree_b),
                "mean_endpoint_degree_centrality": (
                    mean_degree / (node_count - 1) if node_count > 1 else 0
                ),
                "mean_log1p_endpoint_degree": mean_log_degree,
            }
        )

    sl_values = degree_by_label[1]
    non_sl_values = degree_by_label[0]
    if not sl_values or not non_sl_values:
        raise ValueError("Both SL and NONSL candidate pairs are required for the test")

    test = mannwhitneyu(
        sl_values, non_sl_values, alternative="two-sided", method="asymptotic"
    )
    pair_combinations = len(sl_values) * len(non_sl_values)
    probability_sl_greater = test.statistic / pair_combinations
    rank_biserial = 2 * probability_sl_greater - 1
    observed_difference = statistics.mean(sl_values) - statistics.mean(non_sl_values)

    ordered_genes = sorted(candidate_genes, key=int)
    degree_vector = np.asarray(
        [degree_by_gene.get(gene_id, 0) for gene_id in ordered_genes],
        dtype=np.float64,
    )
    sl_incidence = np.zeros(len(ordered_genes), dtype=np.float64)
    non_sl_incidence = np.zeros(len(ordered_genes), dtype=np.float64)
    gene_index = {gene_id: index for index, gene_id in enumerate(ordered_genes)}
    for (gene_a, gene_b), label in selected_pairs:
        incidence = sl_incidence if label == 1 else non_sl_incidence
        incidence[gene_index[gene_a]] += 1
        incidence[gene_index[gene_b]] += 1

    rng = np.random.default_rng(random_seed)
    extreme_permutations = 0
    for _ in range(permutation_count):
        permuted_degrees = rng.permutation(degree_vector)
        permuted_sl_mean = (
            0.5 * np.dot(permuted_degrees, sl_incidence) / len(sl_values)
        )
        permuted_non_sl_mean = (
            0.5 * np.dot(permuted_degrees, non_sl_incidence) / len(non_sl_values)
        )
        if abs(permuted_sl_mean - permuted_non_sl_mean) >= abs(
            observed_difference
        ):
            extreme_permutations += 1
    permutation_p_value = (extreme_permutations + 1) / (permutation_count + 1)

    summary = {
        "string_nodes_in_candidate_universe": node_count,
        "SL_ONLY_pairs_tested": len(sl_values),
        "NONSL_ONLY_pairs_tested": len(non_sl_values),
        "SL_mean_pair_endpoint_degree": statistics.mean(sl_values),
        "NONSL_mean_pair_endpoint_degree": statistics.mean(non_sl_values),
        "SL_median_pair_endpoint_degree": statistics.median(sl_values),
        "NONSL_median_pair_endpoint_degree": statistics.median(non_sl_values),
        "mann_whitney_u": float(test.statistic),
        "mann_whitney_two_sided_p_value": float(test.pvalue),
        "observed_SL_minus_NONSL_mean_degree": observed_difference,
        "node_label_permutation_count": permutation_count,
        "node_label_permutation_seed": random_seed,
        "node_label_permutation_two_sided_p_value": permutation_p_value,
        "probability_random_SL_pair_has_higher_mean_degree": (
            probability_sl_greater
        ),
        "rank_biserial_effect_size": rank_biserial,
        "interpretation_note": (
            "The node-label permutation preserves repeated gene use across pairs "
            "but assumes STRING degree labels are exchangeable across candidate "
            "genes. The Mann-Whitney p-value treats pairs as independent."
        ),
    }
    return pair_rows, summary


def write_plots(
    pair_rows: list[dict], degree_by_gene: dict[str, int], universe_count: int
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sl_values = [
        row["mean_endpoint_degree"] for row in pair_rows if row["label"] == 1
    ]
    non_sl_values = [
        row["mean_endpoint_degree"] for row in pair_rows if row["label"] == 0
    ]
    fig, axis = plt.subplots(figsize=(8, 5))
    axis.boxplot(
        [sl_values, non_sl_values],
        tick_labels=["SL_ONLY", "NONSL_ONLY"],
        showfliers=False,
    )
    axis.set_ylabel("Mean STRING degree of pair endpoints")
    axis.set_title(
        f"Pair endpoint degree in STRING induced network (n={universe_count:,})"
    )
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "string_degree_by_pair_label.png", dpi=180)
    plt.close(fig)

    values = sorted(degree_by_gene.values())
    fig, axis = plt.subplots(figsize=(8, 5))
    positive = [value for value in values if value > 0]
    axis.hist(positive, bins=50)
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlabel("STRING degree (log scale)")
    axis.set_ylabel("Candidate genes (log scale)")
    axis.set_title("STRING degree in feature-supported gene universe")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "string_gene_degree_distribution.png", dpi=180)
    plt.close(fig)


def main() -> int:
    args = parse_args()
    if not 0 <= args.string_min_score <= 1000:
        raise ValueError("--string-min-score must be between 0 and 1000")
    if args.permutations < 1:
        raise ValueError("--permutations must be at least 1")

    required = [
        UNIVERSE_FILE,
        PAIR_SUMMARY_FILE,
        BASELINE_FILE,
        RAW_SL_FILE,
        RAW_NONSL_FILE,
        STRING_LINKS_FILE,
        STRING_ALIASES_FILE,
    ]
    require_files(required)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    gene_symbols, candidate_genes = read_gene_universe(UNIVERSE_FILE)
    raw_sl_rows, raw_sl_pairs = read_raw_pair_file(RAW_SL_FILE)
    raw_nonsl_rows, raw_nonsl_pairs = read_raw_pair_file(RAW_NONSL_FILE)
    all_labels, all_label_counts = read_labeled_pairs(PAIR_SUMMARY_FILE)
    baseline_pairs = read_baseline_pairs(BASELINE_FILE)

    raw_both = raw_sl_pairs & raw_nonsl_pairs
    candidate_baseline_pairs = {
        pair: label
        for pair, label in baseline_pairs.items()
        if pair[0] in candidate_genes and pair[1] in candidate_genes
    }
    candidate_all_label_counts = Counter()
    for pair, label in all_labels.items():
        if pair[0] in candidate_genes and pair[1] in candidate_genes:
            candidate_all_label_counts[
                "SL_ONLY" if label == 1 else "NONSL_ONLY"
            ] += 1
    for row in read_tsv(PAIR_SUMMARY_FILE):
        pair = canonical_pair(row.get("gene_a_id", ""), row.get("gene_b_id", ""))
        if pair and pair[0] in candidate_genes and pair[1] in candidate_genes:
            if row.get("label_status") == "BOTH_LABELS":
                candidate_all_label_counts["BOTH_LABELS"] += 1

    possible_candidate_pairs = len(candidate_genes) * (len(candidate_genes) - 1) // 2
    labeled_candidate_pairs = sum(candidate_all_label_counts.values())

    raw_partners = partner_counts(raw_sl_pairs)
    baseline_sl_pairs = {
        pair for pair, label in baseline_pairs.items() if label == 1
    }
    baseline_sl_partners = partner_counts(baseline_sl_pairs)
    candidate_sl_pairs = {
        pair
        for pair, label in candidate_baseline_pairs.items()
        if label == 1
    }
    candidate_sl_partners = partner_counts(candidate_sl_pairs)

    distribution_rows = [
        {
            "dataset_scope": "raw_synlethdb_union",
            "pair_status": "SL",
            "unique_pairs": len(raw_sl_pairs),
        },
        {
            "dataset_scope": "raw_synlethdb_union",
            "pair_status": "NONSL",
            "unique_pairs": len(raw_nonsl_pairs),
        },
        {
            "dataset_scope": "raw_synlethdb_union",
            "pair_status": "BOTH_LABELS",
            "unique_pairs": len(raw_both),
        },
        *[
            {
                "dataset_scope": "baseline_all_genes",
                "pair_status": status,
                "unique_pairs": all_label_counts[status],
            }
            for status in ("SL_ONLY", "NONSL_ONLY", "BOTH_LABELS")
        ],
        *[
            {
                "dataset_scope": "candidate_feature_universe",
                "pair_status": status,
                "unique_pairs": candidate_all_label_counts[status],
            }
            for status in ("SL_ONLY", "NONSL_ONLY", "BOTH_LABELS")
        ],
    ]
    write_tsv(
        OUTPUT_DIR / "sl_pair_distribution.tsv",
        ["dataset_scope", "pair_status", "unique_pairs"],
        distribution_rows,
    )

    summary_rows = [
        degree_summary(
            "raw_SL_observed_gene_network",
            {gene: len(partners) for gene, partners in raw_partners.items()},
        ),
        degree_summary(
            "baseline_SL_ONLY_all_observed_genes",
            {
                gene: len(partners)
                for gene, partners in baseline_sl_partners.items()
            },
        ),
        degree_summary(
            "baseline_SL_ONLY_candidate_universe",
            {
                gene: len(partners)
                for gene, partners in candidate_sl_partners.items()
            },
            universe_size=len(candidate_genes),
        ),
    ]
    write_tsv(
        OUTPUT_DIR / "sl_partner_distribution_summary.tsv",
        list(summary_rows[0].keys()),
        summary_rows,
    )

    distribution_bin_rows = []
    distribution_bin_rows.extend(
        distribution_bins(
            "raw_SL_observed_gene_network",
            {gene: len(partners) for gene, partners in raw_partners.items()},
        )
    )
    distribution_bin_rows.extend(
        distribution_bins(
            "baseline_SL_ONLY_candidate_universe",
            {
                gene: len(partners)
                for gene, partners in candidate_sl_partners.items()
            },
            universe_size=len(candidate_genes),
        )
    )
    write_tsv(
        OUTPUT_DIR / "sl_partner_count_bins.tsv",
        ["scope", "partner_count_bin", "gene_count"],
        distribution_bin_rows,
    )

    all_genes_for_partner_table = (
        set(gene_symbols) | set(raw_partners) | set(baseline_sl_partners)
    )
    partner_rows = []
    for gene_id in sorted(all_genes_for_partner_table, key=int):
        partner_rows.append(
            {
                "entrez_id": gene_id,
                "gene_symbol": gene_symbols.get(gene_id, ""),
                "in_candidate_universe": gene_id in candidate_genes,
                "raw_SL_unique_partners": len(raw_partners.get(gene_id, set())),
                "baseline_SL_ONLY_unique_partners": len(
                    baseline_sl_partners.get(gene_id, set())
                ),
                "candidate_universe_SL_ONLY_partners": len(
                    candidate_sl_partners.get(gene_id, set())
                ),
            }
        )
    partner_columns = list(partner_rows[0].keys())
    write_tsv(
        OUTPUT_DIR / "sl_partner_counts_all_genes.tsv",
        partner_columns,
        partner_rows,
    )
    for column, filename in [
        ("raw_SL_unique_partners", "top_sl_partners_raw.tsv"),
        (
            "baseline_SL_ONLY_unique_partners",
            "top_sl_partners_baseline.tsv",
        ),
        (
            "candidate_universe_SL_ONLY_partners",
            "top_sl_partners_candidate_universe.tsv",
        ),
    ]:
        top_rows = sorted(
            partner_rows,
            key=lambda row: (-row[column], row["gene_symbol"], int(row["entrez_id"])),
        )[:100]
        write_tsv(OUTPUT_DIR / filename, partner_columns, top_rows)

    write_distribution_plot(
        Counter(
            {
                "SL_ONLY": all_label_counts["SL_ONLY"],
                "NONSL_ONLY": all_label_counts["NONSL_ONLY"],
                "BOTH_LABELS": all_label_counts["BOTH_LABELS"],
            }
        ),
        {
            "raw_sl": {
                gene: len(partners) for gene, partners in raw_partners.items()
            },
            "candidate_sl": {
                gene: len(partners)
                for gene, partners in candidate_sl_partners.items()
            },
        },
        len(candidate_genes),
    )

    query_genes = {
        gene_id
        for pair in candidate_sl_pairs
        for gene_id in pair
    }
    write_enrichment_query(query_genes, candidate_genes)
    enrichment_count = 0
    enrichment_error = ""
    if args.skip_enrichment:
        enrichment_error = "Skipped by --skip-enrichment"
    else:
        enrichment_count, enrichment_error = run_gprofiler(
            query_genes, candidate_genes
        )

    string_degrees, string_metadata = build_string_gene_degrees(
        candidate_genes, args.string_min_score
    )
    centrality_rows, centrality_summary = centrality_association(
        candidate_baseline_pairs,
        candidate_genes,
        string_degrees,
        args.permutations,
    )
    write_tsv(
        OUTPUT_DIR / "string_gene_degree_centrality.tsv",
        [
            "entrez_id",
            "gene_symbol",
            "string_degree",
            "degree_centrality",
        ],
        (
            {
                "entrez_id": gene_id,
                "gene_symbol": gene_symbols.get(gene_id, ""),
                "string_degree": string_degrees[gene_id],
                "degree_centrality": (
                    string_degrees[gene_id] / (len(candidate_genes) - 1)
                    if len(candidate_genes) > 1
                    else 0
                ),
            }
            for gene_id in sorted(candidate_genes, key=int)
        ),
    )
    write_tsv(
        OUTPUT_DIR / "string_degree_by_labeled_pair.tsv",
        list(centrality_rows[0].keys()),
        centrality_rows,
    )
    write_tsv(
        OUTPUT_DIR / "string_centrality_association.tsv",
        ["metric", "value"],
        (
            {"metric": key, "value": value}
            for key, value in {
                **string_metadata,
                **centrality_summary,
            }.items()
        ),
    )
    write_plots(centrality_rows, string_degrees, len(candidate_genes))

    report_lines = [
        "Synthetic lethality dataset analysis",
        "====================================",
        "",
        "Gene universe",
        f"Candidate genes with both DepMap feature types: {len(candidate_genes):,}",
        "This is not every annotated Homo sapiens gene; it is the feature-supported universe.",
        f"All unordered non-self candidate pairs: {possible_candidate_pairs:,}",
        f"Candidate pairs with any SynLethDB label: {labeled_candidate_pairs:,}",
        (
            "Candidate pairs without an observed label: "
            f"{possible_candidate_pairs - labeled_candidate_pairs:,} "
            "(unknown, not confirmed NONSL)"
        ),
        "",
        "Pair-label distributions",
        f"Raw SL rows: {raw_sl_rows:,}; unique unordered pairs: {len(raw_sl_pairs):,}",
        f"Raw NONSL rows: {raw_nonsl_rows:,}; unique unordered pairs: {len(raw_nonsl_pairs):,}",
        f"Pairs in both raw tables: {len(raw_both):,}",
        (
            "Baseline all genes: "
            f"{all_label_counts['SL_ONLY']:,} SL_ONLY, "
            f"{all_label_counts['NONSL_ONLY']:,} NONSL_ONLY, "
            f"{all_label_counts['BOTH_LABELS']:,} BOTH_LABELS"
        ),
        (
            "Baseline pairs within candidate universe: "
            f"{candidate_all_label_counts['SL_ONLY']:,} SL_ONLY, "
            f"{candidate_all_label_counts['NONSL_ONLY']:,} NONSL_ONLY, "
            f"{candidate_all_label_counts['BOTH_LABELS']:,} BOTH_LABELS"
        ),
        "",
        "SL partner counts",
    ]
    for row in summary_rows:
        report_lines.append(
            f"{row['scope']}: genes={row['gene_count_in_summary']:,}, "
            f"with_partner={row['genes_with_at_least_one_partner']:,}, "
            f"mean={row['mean_partners']:.3f}, "
            f"median={row['median_partners']}, "
            f"p90={row['p90_nearest_rank']}, "
            f"p95={row['p95_nearest_rank']}, "
            f"max={row['maximum_partners']}"
        )
    report_lines.extend(
        [
            "",
            "GO/pathway enrichment",
            (
                f"Query: {len(query_genes):,} unique genes that are endpoints "
                "of candidate-universe SL_ONLY pairs"
            ),
            f"Custom background: {len(candidate_genes):,} candidate genes",
            f"g:Profiler results: {enrichment_count:,}",
            (
                f"Enrichment status: ERROR/SKIPPED - {enrichment_error}"
                if enrichment_error
                else "Enrichment status: successful; see go_pathway_enrichment.tsv"
            ),
        ]
    )
    if enrichment_count:
        with (OUTPUT_DIR / "go_pathway_enrichment.tsv").open(
                "r", encoding="utf-8-sig", newline=""
        ) as file:
                enrichment_rows = list(csv.DictReader(file, delimiter="\t"))
        report_lines.extend(["", "Selected significant enrichment terms"])
        excluded_generic_terms = {
                "binding",
                "molecular_function",
                "biological_process",
                "cellular_component",
        }
        for source in GPROFILER_SOURCES:
                terms = [
                    row
                    for row in enrichment_rows
                    if row["source"] == source
                    and "root" not in row["name"].lower()
                    and row["name"].lower() not in excluded_generic_terms
                ]
                terms.sort(key=lambda row: float(row["p_value_fdr"]))
                report_lines.append(f"{source}:")
                for term in terms[:5]:
                    report_lines.append(
                        f"  {term['name']} "
                        f"(FDR p={float(term['p_value_fdr']):.3g}, "
                        f"overlap={term['intersection_size']})"
                    )
    report_lines.extend(
        [
                "",
                "STRING centrality association",
            f"STRING threshold: combined_score >= {args.string_min_score}",
            (
                "STRING graph was mapped to Entrez IDs and induced on the "
                "candidate gene universe."
            ),
            (
                f"SL pairs tested: {centrality_summary['SL_ONLY_pairs_tested']:,}; "
                f"NONSL pairs tested: {centrality_summary['NONSL_ONLY_pairs_tested']:,}"
            ),
            (
                "Mean endpoint degree (SL vs NONSL): "
                f"{centrality_summary['SL_mean_pair_endpoint_degree']:.3f} vs "
                f"{centrality_summary['NONSL_mean_pair_endpoint_degree']:.3f}"
            ),
            (
                "Mann-Whitney pair-independence p-value (naive): "
                f"{centrality_summary['mann_whitney_two_sided_p_value']:.6g}"
            ),
            (
                "Node-label permutation p-value "
                f"({args.permutations:,} permutations): "
                f"{centrality_summary['node_label_permutation_two_sided_p_value']:.6g}"
            ),
            (
                "Rank-biserial effect size: "
                f"{centrality_summary['rank_biserial_effect_size']:.6f}"
            ),
            (
                "The node-label permutation preserves repeated gene use across "
                "pairs but assumes degree labels are exchangeable across genes."
            ),
            (
                "This association is not causal and may reflect study/popularity "
                "bias in the label database."
            ),
        ]
    )
    (OUTPUT_DIR / "analysis_report.txt").write_text(
        "\n".join(report_lines) + "\n", encoding="utf-8"
    )

    print("\n".join(report_lines))
    print(f"\nOutputs written to: {OUTPUT_DIR}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"Analysis failed: {type(error).__name__}: {error}", file=sys.stderr)
        raise