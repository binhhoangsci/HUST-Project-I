from collections import defaultdict
from pathlib import Path
import csv


DATA_DIR = Path("Data")
OUTPUT_DIR = DATA_DIR / "processed_labels"

INPUTS = [
    ("SL", 1, DATA_DIR / "Human.SL.detailed.tsv"),
    ("NONSL", 0, DATA_DIR / "Human.non.SL.detailed.tsv"),
]


def read_tsv(path):
    with path.open("r", encoding="utf-8-sig", newline="") as file:
        yield from csv.DictReader(file, delimiter="\t")


def joined(values):
    return " || ".join(sorted(value for value in values if value))


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    pairs = defaultdict(
        lambda: {
            "labels": set(),
            "sl_rows": 0,
            "nonsl_rows": 0,
            "gene_a_names": set(),
            "gene_b_names": set(),
            "sl_cell_lines": set(),
            "nonsl_cell_lines": set(),
            "sl_sources": set(),
            "nonsl_sources": set(),
            "sl_pmids": set(),
            "nonsl_pmids": set(),
        }
    )

    input_row_counts = {"SL": 0, "NONSL": 0}
    skipped_rows = {"missing_gene_id": 0, "self_pair": 0}

    for table_name, label, path in INPUTS:
        if not path.exists():
            raise FileNotFoundError(f"Không tìm thấy file: {path}")

        for row in read_tsv(path):
            input_row_counts[table_name] += 1

            gene_x_id = (row.get("x:START_ID") or "").strip()
            gene_y_id = (row.get("y:END_ID") or "").strip()
            gene_x_name = (row.get("x_name") or "").strip()
            gene_y_name = (row.get("y_name") or "").strip()

            if not gene_x_id or not gene_y_id:
                skipped_rows["missing_gene_id"] += 1
                continue

            if gene_x_id == gene_y_id:
                skipped_rows["self_pair"] += 1
                continue

            # Mỗi cặp luôn có một thứ tự ID cố định, bất kể đầu vào là A-B hay B-A.
            if gene_x_id <= gene_y_id:
                gene_a_id, gene_b_id = gene_x_id, gene_y_id
                gene_a_name, gene_b_name = gene_x_name, gene_y_name
            else:
                gene_a_id, gene_b_id = gene_y_id, gene_x_id
                gene_a_name, gene_b_name = gene_y_name, gene_x_name

            pair = (gene_a_id, gene_b_id)
            item = pairs[pair]
            item["labels"].add(label)

            if gene_a_name:
                item["gene_a_names"].add(gene_a_name)
            if gene_b_name:
                item["gene_b_names"].add(gene_b_name)

            cell_line = (row.get("cell_line") or "").strip()
            source = (row.get("rel_source") or "").strip()
            pmid = (row.get("pubmed_id") or "").strip()

            if table_name == "SL":
                item["sl_rows"] += 1
                if cell_line:
                    item["sl_cell_lines"].add(cell_line)
                if source:
                    item["sl_sources"].add(source)
                if pmid:
                    item["sl_pmids"].add(pmid)
            else:
                item["nonsl_rows"] += 1
                if cell_line:
                    item["nonsl_cell_lines"].add(cell_line)
                if source:
                    item["nonsl_sources"].add(source)
                if pmid:
                    item["nonsl_pmids"].add(pmid)

    columns = [
        "gene_a_id",
        "gene_a_name",
        "gene_b_id",
        "gene_b_name",
        "label_status",
        "label",
        "sl_evidence_rows",
        "nonsl_evidence_rows",
        "sl_cell_lines",
        "nonsl_cell_lines",
        "sl_sources",
        "nonsl_sources",
        "sl_pubmed_ids",
        "nonsl_pubmed_ids",
    ]

    all_pairs_path = OUTPUT_DIR / "pair_summary_all.tsv"
    baseline_path = OUTPUT_DIR / "pair_baseline_no_conflicts.tsv"
    conflict_path = OUTPUT_DIR / "pairs_with_both_labels.tsv"

    counts = {"SL_ONLY": 0, "NONSL_ONLY": 0, "BOTH_LABELS": 0}
    baseline_rows = 0

    with (
        all_pairs_path.open("w", encoding="utf-8", newline="") as all_file,
        baseline_path.open("w", encoding="utf-8", newline="") as baseline_file,
        conflict_path.open("w", encoding="utf-8", newline="") as conflict_file,
    ):
        all_writer = csv.DictWriter(all_file, fieldnames=columns, delimiter="\t")
        baseline_writer = csv.DictWriter(
            baseline_file, fieldnames=columns, delimiter="\t"
        )
        conflict_writer = csv.DictWriter(
            conflict_file, fieldnames=columns, delimiter="\t"
        )

        all_writer.writeheader()
        baseline_writer.writeheader()
        conflict_writer.writeheader()

        for (gene_a_id, gene_b_id), item in sorted(pairs.items()):
            if item["labels"] == {1}:
                status, label = "SL_ONLY", 1
            elif item["labels"] == {0}:
                status, label = "NONSL_ONLY", 0
            else:
                status, label = "BOTH_LABELS", ""

            counts[status] += 1

            result = {
                "gene_a_id": gene_a_id,
                "gene_a_name": joined(item["gene_a_names"]),
                "gene_b_id": gene_b_id,
                "gene_b_name": joined(item["gene_b_names"]),
                "label_status": status,
                "label": label,
                "sl_evidence_rows": item["sl_rows"],
                "nonsl_evidence_rows": item["nonsl_rows"],
                "sl_cell_lines": joined(item["sl_cell_lines"]),
                "nonsl_cell_lines": joined(item["nonsl_cell_lines"]),
                "sl_sources": joined(item["sl_sources"]),
                "nonsl_sources": joined(item["nonsl_sources"]),
                "sl_pubmed_ids": joined(item["sl_pmids"]),
                "nonsl_pubmed_ids": joined(item["nonsl_pmids"]),
            }

            all_writer.writerow(result)

            if status == "BOTH_LABELS":
                conflict_writer.writerow(result)
            else:
                baseline_writer.writerow(result)
                baseline_rows += 1

    print("Số dòng đầu vào:", input_row_counts)
    print("Số dòng bị bỏ:", skipped_rows)
    print("Số cặp duy nhất:", len(pairs))
    print("Cặp chỉ có nhãn SL:", counts["SL_ONLY"])
    print("Cặp chỉ có nhãn NONSL:", counts["NONSL_ONLY"])
    print("Cặp có cả hai nhãn:", counts["BOTH_LABELS"])
    print("Số dòng trong baseline:", baseline_rows)
    print("Đã ghi kết quả vào:", OUTPUT_DIR)


if __name__ == "__main__":
    main()