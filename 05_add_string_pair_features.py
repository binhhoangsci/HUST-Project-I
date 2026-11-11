import csv
import gzip
import json
import sqlite3
import tempfile
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "Data"
INPUT_DIR = DATA_DIR / "model_dataset"
OUTPUT_DIR = DATA_DIR / "model_dataset_final"
GENE_SUMMARY_FILE = INPUT_DIR / "gene_level_depmap_summary.tsv"
STRING_LINKS_FILE = DATA_DIR / "9606.protein.links.v12.5.txt.gz"
STRING_ALIASES_FILE = DATA_DIR / "9606.protein.aliases.v12.5.txt.gz"
STRING_GENE_ID_SOURCE = "Ensembl_HGNC_entrez_id"
SPLITS = ("train", "validation", "test")
STRING_FEATURE_COLUMNS = ["string_score", "string_has_recorded_edge"]


def read_tsv(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as file:
        yield from csv.DictReader(file, delimiter="\t")


def write_tsv(path: Path, columns: list[str], rows) -> None:
    with path.open("w", encoding="utf-8", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def ordered_pair(gene_a: str, gene_b: str) -> tuple[str, str]:
    return tuple(sorted((gene_a, gene_b), key=int))


def read_candidate_genes(path: Path) -> tuple[set[str], dict[str, str]]:
    gene_ids = set()
    gene_names = {}
    for row in read_tsv(path):
        gene_id = row["entrez_id"].strip()
        if not gene_id:
            raise ValueError(f"Blank Entrez ID in {path}")
        gene_ids.add(gene_id)
        gene_names[gene_id] = row["gene_name"].strip() or gene_id
    return gene_ids, gene_names


def map_string_proteins(
    candidate_genes: set[str],
) -> tuple[dict[str, str], dict[str, int], int]:
    protein_gene_candidates: dict[str, set[str]] = defaultdict(set)
    with gzip.open(STRING_ALIASES_FILE, "rt", encoding="utf-8") as file:
        reader = csv.DictReader(file, delimiter="\t")
        required = {"#string_protein_id", "alias", "source"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"STRING aliases file is missing columns: {sorted(missing)}"
            )
        for row in reader:
            if row["source"] != STRING_GENE_ID_SOURCE:
                continue
            gene_id = row["alias"].strip()
            protein_id = row["#string_protein_id"].strip()
            if gene_id in candidate_genes and protein_id:
                protein_gene_candidates[protein_id].add(gene_id)

    protein_to_gene = {
        protein_id: next(iter(gene_ids))
        for protein_id, gene_ids in protein_gene_candidates.items()
        if len(gene_ids) == 1
    }
    ambiguous_protein_count = sum(
        len(gene_ids) > 1 for gene_ids in protein_gene_candidates.values()
    )
    proteins_per_gene: Counter = Counter(protein_to_gene.values())
    return protein_to_gene, dict(proteins_per_gene), ambiguous_protein_count


def build_gene_pair_score_db(
    candidate_genes: set[str],
    protein_to_gene: dict[str, str],
    database_path: Path,
) -> dict[str, int]:
    metrics: Counter = Counter()
    connection = sqlite3.connect(database_path)
    try:
        connection.execute(
            "CREATE TABLE edge ("
            "gene_a TEXT NOT NULL, "
            "gene_b TEXT NOT NULL, "
            "score INTEGER NOT NULL, "
            "PRIMARY KEY (gene_a, gene_b)"
            ") WITHOUT ROWID"
        )
        insert_sql = (
            "INSERT INTO edge (gene_a, gene_b, score) VALUES (?, ?, ?) "
            "ON CONFLICT(gene_a, gene_b) DO UPDATE "
            "SET score = MAX(edge.score, excluded.score)"
        )
        batch = []

        with gzip.open(STRING_LINKS_FILE, "rt", encoding="utf-8") as file:
            header = next(file, "").split()
            if header != ["protein1", "protein2", "combined_score"]:
                raise ValueError(
                    f"Unexpected STRING links header: {header}"
                )
            for line_number, line in enumerate(file, start=2):
                metrics["input_protein_links"] += 1
                fields = line.split()
                if len(fields) != 3:
                    raise ValueError(
                        f"Malformed STRING links row {line_number}"
                    )
                protein_a, protein_b, score_text = fields
                try:
                    score = int(score_text)
                except ValueError as error:
                    raise ValueError(
                        f"Invalid STRING combined_score at row {line_number}"
                    ) from error
                if not 0 <= score <= 1000:
                    raise ValueError(
                        f"STRING score outside 0-1000 at row {line_number}"
                    )

                gene_a = protein_to_gene.get(protein_a)
                gene_b = protein_to_gene.get(protein_b)
                if gene_a is None or gene_b is None:
                    metrics["links_with_unmapped_protein"] += 1
                    continue
                if gene_a not in candidate_genes or gene_b not in candidate_genes:
                    metrics["links_outside_candidate_universe"] += 1
                    continue
                if gene_a == gene_b:
                    metrics["within_gene_links_excluded"] += 1
                    continue

                pair = ordered_pair(gene_a, gene_b)
                batch.append((pair[0], pair[1], score))
                metrics["mapped_protein_links"] += 1
                if len(batch) >= 50000:
                    connection.executemany(insert_sql, batch)
                    connection.commit()
                    batch.clear()

        if batch:
            connection.executemany(insert_sql, batch)
            connection.commit()

        metrics["unique_gene_pairs_after_mapping"] = connection.execute(
            "SELECT COUNT(*) FROM edge"
        ).fetchone()[0]
        connection.execute("CREATE INDEX edge_score_idx ON edge(score)")
        connection.commit()
        return dict(metrics)
    finally:
        connection.close()


def read_split_rows(path: Path) -> tuple[list[dict], list[str]]:
    rows = list(read_tsv(path))
    if not rows:
        raise ValueError(f"Empty split file: {path}")
    return rows, list(rows[0].keys())


def main() -> None:
    required = [
        GENE_SUMMARY_FILE,
        STRING_LINKS_FILE,
        STRING_ALIASES_FILE,
        *(INPUT_DIR / f"{split}.tsv" for split in SPLITS),
        INPUT_DIR / "pair_split_assignments.tsv",
    ]
    for path in required:
        if not path.is_file():
            raise FileNotFoundError(f"Required input not found: {path}")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    candidate_genes, gene_names = read_candidate_genes(GENE_SUMMARY_FILE)
    protein_to_gene, proteins_per_gene, ambiguous_proteins = map_string_proteins(
        candidate_genes
    )

    with tempfile.TemporaryDirectory(
        prefix="string_pair_features_", dir=OUTPUT_DIR
    ) as temporary_directory:
        database_path = Path(temporary_directory) / "string_scores.sqlite"
        mapping_metrics = build_gene_pair_score_db(
            candidate_genes, protein_to_gene, database_path
        )

        mapping_rows = [
            {
                "entrez_id": gene_id,
                "gene_symbol": gene_names[gene_id],
                "mapped_string_protein_count": proteins_per_gene.get(gene_id, 0),
                "has_string_mapping": int(proteins_per_gene.get(gene_id, 0) > 0),
            }
            for gene_id in sorted(candidate_genes, key=int)
        ]
        write_tsv(
            OUTPUT_DIR / "string_gene_mapping_audit.tsv",
            [
                "entrez_id",
                "gene_symbol",
                "mapped_string_protein_count",
                "has_string_mapping",
            ],
            mapping_rows,
        )

        connection = sqlite3.connect(database_path)
        try:
            pair_score_query = (
                "SELECT score FROM edge WHERE gene_a = ? AND gene_b = ?"
            )
            global_pair_ids = set()
            split_label_counts = {}
            split_rows_written = {}
            output_columns = None

            for split in SPLITS:
                input_rows, input_columns = read_split_rows(
                    INPUT_DIR / f"{split}.tsv"
                )
                if output_columns is None:
                    output_columns = input_columns + STRING_FEATURE_COLUMNS
                elif input_columns != output_columns[:-len(STRING_FEATURE_COLUMNS)]:
                    raise ValueError(
                        f"Input columns differ across split {split}"
                    )

                output_rows = []
                labels = Counter()
                local_pair_ids = set()
                recorded_edges = 0
                for row in input_rows:
                    gene_a = row["gene_a_id"].strip()
                    gene_b = row["gene_b_id"].strip()
                    pair = ordered_pair(gene_a, gene_b)
                    expected_pair_id = f"{pair[0]}_{pair[1]}"
                    if row["pair_id"] != expected_pair_id:
                        raise ValueError(
                            f"Pair ID mismatch in {split}: {row['pair_id']}"
                        )
                    if pair[0] not in candidate_genes or pair[1] not in candidate_genes:
                        raise ValueError(
                            f"Pair lacks candidate feature gene in {split}: {pair}"
                        )
                    if row["label"] not in {"0", "1"}:
                        raise ValueError(
                            f"Invalid label in {split}: {row['label']}"
                        )
                    if row["label_status"] not in {"SL_ONLY", "NONSL_ONLY"}:
                        raise ValueError(
                            f"Unexpected label status in {split}: "
                            f"{row['label_status']}"
                        )
                    if pair in local_pair_ids:
                        raise ValueError(f"Duplicate pair in {split}: {pair}")
                    local_pair_ids.add(pair)

                    result = connection.execute(
                        pair_score_query, pair
                    ).fetchone()
                    has_recorded_edge = int(result is not None)
                    score = result[0] / 1000 if result is not None else 0.0
                    recorded_edges += has_recorded_edge
                    labels[int(row["label"])] += 1
                    output_rows.append(
                        {
                            **row,
                            "string_score": score,
                            "string_has_recorded_edge": has_recorded_edge,
                        }
                    )

                overlap = global_pair_ids & local_pair_ids
                if overlap:
                    raise ValueError(
                        f"Pairs duplicated across splits: {list(overlap)[:5]}"
                    )
                global_pair_ids.update(local_pair_ids)
                write_tsv(
                    OUTPUT_DIR / f"{split}.tsv",
                    output_columns,
                    output_rows,
                )
                split_label_counts[split] = labels
                split_rows_written[split] = {
                    "rows": len(output_rows),
                    "SL_label_1": labels[1],
                    "NONSL_label_0": labels[0],
                    "pairs_with_STRING_record": recorded_edges,
                    "pairs_without_STRING_record": len(output_rows) - recorded_edges,
                }

            expected_assignment_count = sum(
                1 for _ in read_tsv(INPUT_DIR / "pair_split_assignments.tsv")
            )
            if len(global_pair_ids) != expected_assignment_count:
                raise ValueError(
                    "Combined output pair count does not match the saved "
                    "split assignments"
                )
            write_tsv(
                OUTPUT_DIR / "split_summary.tsv",
                [
                    "split",
                    "rows",
                    "SL_label_1",
                    "NONSL_label_0",
                    "pairs_with_STRING_record",
                    "pairs_without_STRING_record",
                ],
                (
                    {"split": split, **split_rows_written[split]}
                    for split in SPLITS
                ),
            )
        finally:
            connection.close()

    mapping_count = len(protein_to_gene)
    summary = {
        "inputs": {
            "string_links": str(STRING_LINKS_FILE.relative_to(ROOT)),
            "string_aliases": str(STRING_ALIASES_FILE.relative_to(ROOT)),
            "string_version": "12.5",
            "gene_id_alias_source": STRING_GENE_ID_SOURCE,
            "score_range": [0, 1000],
        },
        "gene_mapping": {
            "candidate_genes": len(candidate_genes),
            "candidate_genes_with_string_mapping": len(proteins_per_gene),
            "candidate_genes_without_string_mapping": (
                len(candidate_genes) - len(proteins_per_gene)
            ),
            "unambiguous_string_proteins_mapped": mapping_count,
            "ambiguous_string_proteins_excluded": ambiguous_proteins,
        },
        "link_processing": mapping_metrics,
        "pair_feature_rules": {
            "string_score": (
                "Maximum combined_score/1000 among all STRING protein links "
                "mapping to the unordered gene pair."
            ),
            "string_has_recorded_edge": (
                "1 if at least one STRING link maps to the pair; otherwise 0."
            ),
            "absent_link_value": (
                "Score 0 and edge flag 0 means no link recorded in this STRING "
                "file, not experimentally proven absence of association."
            ),
            "string_edge_type": (
                "STRING combined functional association; not necessarily a "
                "direct physical interaction."
            ),
            "protein_to_gene_collapse": (
                "For multiple protein links mapping to one unordered gene pair, "
                "retain the maximum combined_score."
            ),
        },
        "split_policy": (
            "Reuses existing pair_split_assignments.tsv exactly; no resplitting."
        ),
    }
    (OUTPUT_DIR / "string_preprocessing_metadata.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )

    report_lines = [
        "Final pair-level dataset with STRING features",
        "==============================================",
        f"Candidate genes: {len(candidate_genes):,}",
        (
            "Candidate genes with at least one unambiguous STRING protein "
            f"mapping: {len(proteins_per_gene):,}"
        ),
        f"Unambiguous STRING protein IDs mapped: {mapping_count:,}",
        f"Ambiguous protein mappings excluded: {ambiguous_proteins:,}",
        f"Mapped STRING protein links processed: {mapping_metrics['mapped_protein_links']:,}",
        (
            "Unique gene-pair STRING edges after mapping: "
            f"{mapping_metrics['unique_gene_pairs_after_mapping']:,}"
        ),
        (
            "For each unordered gene pair, STRING score is the maximum "
            "combined_score across mapped protein links, divided by 1000."
        ),
        (
            "No recorded STRING link is encoded as score=0 and "
            "string_has_recorded_edge=0; this is not proof of no biological "
            "association."
        ),
        (
            "STRING combined_score describes functional association and is "
            "not necessarily a direct physical interaction."
        ),
        "Existing train/validation/test assignments were reused unchanged.",
        "",
        "Split counts:",
    ]
    for split, values in split_rows_written.items():
        report_lines.append(
            f"{split}: n={values['rows']:,}, "
            f"SL={values['SL_label_1']:,}, "
            f"NONSL={values['NONSL_label_0']:,}, "
            f"STRING edge={values['pairs_with_STRING_record']:,}, "
            f"no recorded edge={values['pairs_without_STRING_record']:,}"
        )
    (OUTPUT_DIR / "dataset_report.txt").write_text(
        "\n".join(report_lines) + "\n", encoding="utf-8"
    )
    print("\n".join(report_lines))
    print(f"\nFinal dataset written to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
