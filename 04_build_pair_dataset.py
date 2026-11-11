import csv
import math
import random
import statistics
from collections import Counter, defaultdict
from pathlib import Path

from sklearn.model_selection import train_test_split


ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "Data"
LABELS_FILE = DATA_DIR / "processed_labels" / "pair_baseline_no_conflicts.tsv"
CRISPR_FILE = (
    DATA_DIR / "processed_depmap" / "crispr_effect_common_models.tsv"
)
EXPRESSION_FILE = (
    DATA_DIR / "processed_depmap" / "expression_logp1_common_models.tsv"
)
GENE_COVERAGE_FILE = DATA_DIR / "processed_depmap" / "gene_coverage.tsv"
OUTPUT_DIR = DATA_DIR / "model_dataset"
RANDOM_SEED = 42

PAIR_FEATURES = [
    "crispr_min",
    "crispr_max",
    "crispr_mean",
    "crispr_abs_diff",
    "expression_min",
    "expression_max",
    "expression_mean",
    "expression_abs_diff",
]

PAIR_COLUMNS = [
    "pair_id",
    "gene_a_id",
    "gene_a_name",
    "gene_b_id",
    "gene_b_name",
    "label",
    "label_status",
    *PAIR_FEATURES,
]


def read_rows(path: Path, delimiter: str = "\t"):
    with path.open("r", encoding="utf-8-sig", newline="") as file:
        yield from csv.DictReader(file, delimiter=delimiter)


def read_feature_matrix(
    path: Path,
) -> tuple[dict[str, dict[str, float | None]], list[str], dict[str, float]]:
    profiles: dict[str, dict[str, float | None]] = {}
    with path.open("r", encoding="utf-8-sig", newline="") as file:
        reader = csv.reader(file, delimiter="\t")
        header = next(reader)
        if not header or header[0] != "ModelID":
            raise ValueError(f"{path.name}: first column must be ModelID")
        gene_ids = header[1:]
        if len(gene_ids) != len(set(gene_ids)):
            raise ValueError(f"{path.name}: duplicate gene ID columns")

        for line_number, row in enumerate(reader, start=2):
            if len(row) != len(header):
                raise ValueError(
                    f"{path.name}:{line_number}: expected {len(header)} "
                    f"columns, found {len(row)}"
                )
            model_id = row[0].strip()
            if not model_id or model_id in profiles:
                raise ValueError(
                    f"{path.name}:{line_number}: blank or duplicate ModelID "
                    f"{model_id!r}"
                )
            values: dict[str, float | None] = {}
            for gene_id, raw_value in zip(gene_ids, row[1:]):
                raw_value = raw_value.strip()
                if not raw_value:
                    values[gene_id] = None
                    continue
                try:
                    value = float(raw_value)
                except ValueError as error:
                    raise ValueError(
                        f"{path.name}:{line_number}: invalid numeric value "
                        f"for gene {gene_id}: {raw_value!r}"
                    ) from error
                values[gene_id] = value if math.isfinite(value) else None
            profiles[model_id] = values

    return profiles, gene_ids, {}


def summarize_gene_features(
    profiles: dict[str, dict[str, float | None]],
    gene_ids: list[str],
    value_name: str,
) -> dict[str, dict[str, float]]:
    summaries = {}
    profile_count = len(profiles)
    if not profile_count:
        raise ValueError(f"No profiles found for {value_name}")

    for gene_id in gene_ids:
        observed = [
            profile[gene_id]
            for profile in profiles.values()
            if profile[gene_id] is not None
        ]
        if not observed:
            raise ValueError(
                f"Gene {gene_id} has no observed {value_name} values"
            )
        summaries[gene_id] = {
            "median": statistics.median(observed),
            "missing_fraction": (profile_count - len(observed)) / profile_count,
            "observed_profiles": len(observed),
        }
    return summaries


def read_gene_names(path: Path) -> dict[str, str]:
    names: dict[str, set[str]] = defaultdict(set)
    for row in read_rows(path):
        names[row["gene_a_id"]].add(row["gene_a_name"])
        names[row["gene_b_id"]].add(row["gene_b_name"])
    return {
        gene_id: sorted(name for name in gene_names if name)[0]
        if any(gene_names)
        else gene_id
        for gene_id, gene_names in names.items()
    }


def make_symmetric_features(
    gene_a_id: str,
    gene_b_id: str,
    crispr: dict[str, dict[str, float]],
    expression: dict[str, dict[str, float]],
) -> dict[str, float]:
    values = {}
    for prefix, summaries in [
        ("crispr", crispr),
        ("expression", expression),
    ]:
        a = summaries[gene_a_id]["median"]
        b = summaries[gene_b_id]["median"]
        values[f"{prefix}_min"] = min(a, b)
        values[f"{prefix}_max"] = max(a, b)
        values[f"{prefix}_mean"] = (a + b) / 2
        values[f"{prefix}_abs_diff"] = abs(a - b)
    return values


def write_rows(path: Path, columns: list[str], rows) -> None:
    with path.open("w", encoding="utf-8", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    for path in [
        LABELS_FILE,
        CRISPR_FILE,
        EXPRESSION_FILE,
        GENE_COVERAGE_FILE,
    ]:
        if not path.is_file():
            raise FileNotFoundError(f"Required input not found: {path}")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    crispr_profiles, crispr_gene_ids, _ = read_feature_matrix(CRISPR_FILE)
    expression_profiles, expression_gene_ids, _ = read_feature_matrix(
        EXPRESSION_FILE
    )
    if set(crispr_profiles) != set(expression_profiles):
        raise ValueError("CRISPR and expression ModelID sets do not match")
    if set(crispr_gene_ids) != set(expression_gene_ids):
        raise ValueError("CRISPR and expression gene ID sets do not match")

    gene_id_set = set(crispr_gene_ids)
    crispr_summary = summarize_gene_features(
        crispr_profiles, crispr_gene_ids, "CRISPR gene effect"
    )
    expression_summary = summarize_gene_features(
        expression_profiles, expression_gene_ids, "expression"
    )
    gene_names = read_gene_names(LABELS_FILE)

    gene_summary_rows = []
    for gene_id in sorted(gene_id_set, key=int):
        gene_summary_rows.append(
            {
                "entrez_id": gene_id,
                "gene_name": gene_names.get(gene_id, gene_id),
                "crispr_median_across_models": crispr_summary[gene_id]["median"],
                "crispr_missing_fraction": crispr_summary[gene_id][
                    "missing_fraction"
                ],
                "crispr_observed_models": crispr_summary[gene_id][
                    "observed_profiles"
                ],
                "expression_median_across_models": expression_summary[gene_id][
                    "median"
                ],
                "expression_missing_fraction": expression_summary[gene_id][
                    "missing_fraction"
                ],
                "expression_observed_models": expression_summary[gene_id][
                    "observed_profiles"
                ],
            }
        )
    write_rows(
        OUTPUT_DIR / "gene_level_depmap_summary.tsv",
        list(gene_summary_rows[0].keys()),
        gene_summary_rows,
    )

    pair_rows = []
    seen_pairs = set()
    skipped_pairs_missing_features = 0
    for row in read_rows(LABELS_FILE):
        gene_a_id = row["gene_a_id"].strip()
        gene_b_id = row["gene_b_id"].strip()
        pair = tuple(sorted((gene_a_id, gene_b_id), key=int))
        if pair in seen_pairs:
            raise ValueError(f"Duplicate pair in baseline labels: {pair}")
        seen_pairs.add(pair)

        if pair[0] not in gene_id_set or pair[1] not in gene_id_set:
            skipped_pairs_missing_features += 1
            continue
        if row["label"] not in {"0", "1"}:
            raise ValueError(f"Invalid label {row['label']!r} for pair {pair}")

        gene_a_id, gene_b_id = pair
        pair_features = make_symmetric_features(
            gene_a_id, gene_b_id, crispr_summary, expression_summary
        )
        pair_rows.append(
            {
                "pair_id": f"{gene_a_id}_{gene_b_id}",
                "gene_a_id": gene_a_id,
                "gene_a_name": gene_names.get(gene_a_id, gene_a_id),
                "gene_b_id": gene_b_id,
                "gene_b_name": gene_names.get(gene_b_id, gene_b_id),
                "label": int(row["label"]),
                "label_status": row["label_status"],
                **pair_features,
            }
        )

    if not pair_rows:
        raise ValueError("No labeled pairs have features for both genes")

    labels = [row["label"] for row in pair_rows]
    train_rows, temporary_rows = train_test_split(
        pair_rows,
        test_size=0.30,
        random_state=RANDOM_SEED,
        stratify=labels,
    )
    temporary_labels = [row["label"] for row in temporary_rows]
    validation_rows, test_rows = train_test_split(
        temporary_rows,
        test_size=0.50,
        random_state=RANDOM_SEED,
        stratify=temporary_labels,
    )

    partitions = {
        "train": train_rows,
        "validation": validation_rows,
        "test": test_rows,
    }
    all_split_pairs = set()
    split_assignments = []
    split_summary = []
    for split_name, rows in partitions.items():
        pair_keys = {
            (row["gene_a_id"], row["gene_b_id"])
            for row in rows
        }
        if all_split_pairs & pair_keys:
            raise ValueError("The same gene pair occurs in multiple splits")
        all_split_pairs.update(pair_keys)
        write_rows(
            OUTPUT_DIR / f"{split_name}.tsv",
            PAIR_COLUMNS,
            rows,
        )
        counts = Counter(row["label"] for row in rows)
        split_summary.append(
            {
                "split": split_name,
                "rows": len(rows),
                "SL_label_1": counts[1],
                "NONSL_label_0": counts[0],
                "SL_fraction": counts[1] / len(rows),
            }
        )
        split_assignments.extend(
            {"pair_id": row["pair_id"], "split": split_name}
            for row in rows
        )

    write_rows(
        OUTPUT_DIR / "pair_split_summary.tsv",
        ["split", "rows", "SL_label_1", "NONSL_label_0", "SL_fraction"],
        split_summary,
    )
    write_rows(
        OUTPUT_DIR / "pair_split_assignments.tsv",
        ["pair_id", "split"],
        sorted(split_assignments, key=lambda row: row["pair_id"]),
    )

    report_lines = [
        "Pair-level feature dataset",
        "==========================",
        f"Input labeled pairs: {len(seen_pairs)}",
        f"Pairs skipped because one or both genes lack selected features: "
        f"{skipped_pairs_missing_features}",
        f"Pairs retained: {len(pair_rows)}",
        f"DepMap models summarized: {len(crispr_profiles)}",
        f"Genes summarized: {len(gene_id_set)}",
        f"Features per pair: {len(PAIR_FEATURES)}",
        f"Split seed: {RANDOM_SEED}",
        "Split method: stratified random split by unique unordered gene pair",
        "Split proportions: approximately 70% train, 15% validation, 15% test",
        (
            "This is a warm-start split: a gene may appear in more than one "
            "split, but a gene pair is assigned to exactly one split."
        ),
        (
            "CRISPR missing-value strategy: summarize each gene using the "
            "median of observed model values; missing cells are excluded from "
            "the median, not replaced by zero. No gene in this selected matrix "
            "has all model values missing."
        ),
        (
            "CRISPR missing fractions are recorded in gene_level_depmap_summary.tsv "
            "for audit, but are not model input features."
        ),
        (
            "Expression values were already log-transformed in DepMap; no "
            "additional log transform or scaling was applied."
        ),
        (
            "Pair feature construction: for each gene-level median (CRISPR and "
            "expression), compute min, max, mean and absolute difference for "
            "the two genes. These features are invariant to swapping gene order."
        ),
        (
            "No cell-line context is added: this is pair-level prediction, not "
            "pair-plus-cell-line prediction."
        ),
        (
            "No standardization/scaling was fit here. If the selected model "
            "needs scaling, fit it on train only and apply it to validation/test."
        ),
        "",
        "Split counts:",
    ]
    for row in split_summary:
        report_lines.append(
            f"{row['split']}: n={row['rows']}, "
            f"SL={row['SL_label_1']}, NONSL={row['NONSL_label_0']}, "
            f"SL_fraction={row['SL_fraction']:.6f}"
        )
    (OUTPUT_DIR / "dataset_report.txt").write_text(
        "\n".join(report_lines) + "\n", encoding="utf-8"
    )

    print("\n".join(report_lines))
    print(f"\nOutputs written to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
