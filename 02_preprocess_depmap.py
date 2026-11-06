import csv
import math
import re
from collections import defaultdict
from pathlib import Path


DATA_DIR = Path("Data")
DEPMAP_DIR = DATA_DIR / "DepMap"
LABEL_FILE = DATA_DIR / "processed_labels" / "pair_baseline_no_conflicts.tsv"
CRISPR_FILE = DEPMAP_DIR / "CRISPRGeneEffect.csv"
EXPRESSION_FILE = DEPMAP_DIR / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
MODEL_FILE = DEPMAP_DIR / "Model.csv"
OUTPUT_DIR = DATA_DIR / "processed_depmap"

EXPRESSION_METADATA_COLUMNS = 6
MODEL_METADATA_COLUMNS = [
    "ModelID",
    "CellLineName",
    "StrippedCellLineName",
    "DepmapModelType",
    "OncotreeLineage",
    "OncotreePrimaryDisease",
    "OncotreeSubtype",
    "OncotreeCode",
]


def read_header(path, delimiter=","):
    with path.open("r", encoding="utf-8-sig", newline="") as file:
        return next(csv.reader(file, delimiter=delimiter))


def gene_columns(header, start_index, file_label):
    columns = {}
    symbols = {}

    for index in range(start_index, len(header)):
        column_name = header[index]
        match = re.search(r"\((\d+)\)\s*$", column_name)
        if not match:
            raise ValueError(
                f"Không đọc được Entrez ID trong cột {index + 1} "
                f"của {file_label}: {column_name!r}"
            )

        gene_id = match.group(1)
        if gene_id in columns:
            raise ValueError(
                f"Entrez ID {gene_id} bị lặp trong header của {file_label}"
            )

        symbol = column_name[: match.start()].strip()
        columns[gene_id] = index
        symbols[gene_id] = symbol

    return columns, symbols


def read_label_genes(path):
    genes = defaultdict(set)
    with path.open("r", encoding="utf-8-sig", newline="") as file:
        reader = csv.DictReader(file, delimiter="\t")
        required = {
            "gene_a_id",
            "gene_a_name",
            "gene_b_id",
            "gene_b_name",
            "label",
        }
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Thiếu cột trong bảng nhãn: {sorted(missing)}")

        for row in reader:
            if row["label"] not in {"0", "1"}:
                raise ValueError(
                    "Bảng baseline phải chỉ có nhãn 0 hoặc 1; "
                    f"gặp nhãn {row['label']!r}"
                )
            genes[row["gene_a_id"]].add(row["gene_a_name"])
            genes[row["gene_b_id"]].add(row["gene_b_name"])

    return genes


def read_model_metadata(path):
    models = {}
    with path.open("r", encoding="utf-8-sig", newline="") as file:
        reader = csv.DictReader(file)
        missing = set(MODEL_METADATA_COLUMNS).difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Thiếu cột trong Model.csv: {sorted(missing)}")

        for row in reader:
            model_id = row["ModelID"].strip()
            if model_id in models:
                raise ValueError(f"ModelID bị lặp trong Model.csv: {model_id}")
            models[model_id] = row

    return models


def read_crispr_model_ids(path):
    model_ids = set()
    with path.open("r", encoding="utf-8-sig", newline="") as file:
        reader = csv.reader(file)
        next(reader)
        for row in reader:
            if row and row[0].strip():
                model_id = row[0].strip()
                if model_id in model_ids:
                    raise ValueError(
                        f"ModelID bị lặp trong CRISPRGeneEffect.csv: {model_id}"
                    )
                model_ids.add(model_id)
    return model_ids


def read_default_expression_model_ids(path):
    model_ids = set()
    with path.open("r", encoding="utf-8-sig", newline="") as file:
        reader = csv.DictReader(file)
        required = {"ModelID", "IsDefaultEntryForModel"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"Thiếu cột metadata trong file biểu hiện: {sorted(missing)}"
            )

        for row in reader:
            if row["IsDefaultEntryForModel"].strip().lower() == "yes":
                model_id = row["ModelID"].strip()
                if model_id in model_ids:
                    raise ValueError(
                        "Có nhiều hàng expression mặc định cho ModelID "
                        f"{model_id}"
                    )
                model_ids.add(model_id)

    return model_ids


def normalize_numeric(value, path, model_id, gene_id, missing_counts):
    value = value.strip()
    if not value or value.lower() in {"na", "nan", "null"}:
        missing_counts[0] += 1
        return ""

    try:
        number = float(value)
    except ValueError as error:
        raise ValueError(
            f"Giá trị không phải số trong {path}, "
            f"ModelID={model_id}, EntrezID={gene_id}: {value!r}"
        ) from error

    if not math.isfinite(number):
        missing_counts[0] += 1
        return ""

    return value


def write_feature_matrix(
    source_path,
    output_path,
    selected_models,
    selected_gene_ids,
    source_gene_columns,
    model_id_column,
    expression=False,
):
    missing_counts = [0]
    written_models = set()
    selected_indices = [
        source_gene_columns[gene_id] for gene_id in selected_gene_ids
    ]

    with (
        source_path.open("r", encoding="utf-8-sig", newline="") as source_file,
        output_path.open("w", encoding="utf-8", newline="") as output_file,
    ):
        reader = csv.reader(source_file)
        source_header = next(reader)
        writer = csv.writer(output_file, delimiter="\t")
        writer.writerow(["ModelID", *selected_gene_ids])

        for row in reader:
            if not row:
                continue

            model_id = row[model_id_column].strip()
            if model_id not in selected_models:
                continue

            if expression:
                if row[5].strip().lower() != "yes":
                    continue
                if model_id in written_models:
                    raise ValueError(
                        "Có nhiều hàng expression mặc định cho ModelID "
                        f"{model_id}"
                    )

            values = [
                normalize_numeric(
                    row[index],
                    source_path.name,
                    model_id,
                    gene_id,
                    missing_counts,
                )
                for gene_id, index in zip(selected_gene_ids, selected_indices)
            ]
            writer.writerow([model_id, *values])
            written_models.add(model_id)

    if written_models != selected_models:
        missing_models = sorted(selected_models - written_models)
        raise ValueError(
            f"{source_path.name}: thiếu {len(missing_models)} ModelID "
            f"được chọn; ví dụ {missing_models[:5]}"
        )

    return len(written_models), missing_counts[0]


def main():
    required_files = [
        LABEL_FILE,
        CRISPR_FILE,
        EXPRESSION_FILE,
        MODEL_FILE,
    ]
    for path in required_files:
        if not path.exists():
            raise FileNotFoundError(f"Không tìm thấy file: {path}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    label_genes = read_label_genes(LABEL_FILE)
    crispr_header = read_header(CRISPR_FILE)
    expression_header = read_header(EXPRESSION_FILE)
    crispr_columns, crispr_symbols = gene_columns(
        crispr_header, 1, CRISPR_FILE.name
    )
    expression_columns, expression_symbols = gene_columns(
        expression_header, EXPRESSION_METADATA_COLUMNS, EXPRESSION_FILE.name
    )

    model_metadata = read_model_metadata(MODEL_FILE)
    crispr_models = read_crispr_model_ids(CRISPR_FILE)
    expression_models = read_default_expression_model_ids(EXPRESSION_FILE)

    common_models = (
        crispr_models
        & expression_models
        & set(model_metadata)
    )
    if not common_models:
        raise ValueError("Không tìm thấy ModelID chung giữa ba file DepMap")

    selected_gene_ids = sorted(
        set(label_genes) & set(crispr_columns) & set(expression_columns),
        key=lambda gene_id: int(gene_id),
    )
    if not selected_gene_ids:
        raise ValueError(
            "Không tìm thấy gene chung giữa bảng nhãn và hai ma trận DepMap"
        )

    gene_map_path = OUTPUT_DIR / "gene_coverage.tsv"
    with gene_map_path.open("w", encoding="utf-8", newline="") as file:
        writer = csv.writer(file, delimiter="\t")
        writer.writerow(
            [
                "entrez_id",
                "label_gene_names",
                "crispr_symbol",
                "expression_symbol",
                "in_crispr",
                "in_expression",
                "in_both_feature_matrices",
            ]
        )
        for gene_id in sorted(label_genes, key=lambda value: int(value)):
            in_crispr = gene_id in crispr_columns
            in_expression = gene_id in expression_columns
            writer.writerow(
                [
                    gene_id,
                    " || ".join(sorted(name for name in label_genes[gene_id] if name)),
                    crispr_symbols.get(gene_id, ""),
                    expression_symbols.get(gene_id, ""),
                    in_crispr,
                    in_expression,
                    in_crispr and in_expression,
                ]
            )

    model_metadata_path = OUTPUT_DIR / "model_metadata_common.tsv"
    with model_metadata_path.open("w", encoding="utf-8", newline="") as file:
        writer = csv.DictWriter(
            file,
            fieldnames=MODEL_METADATA_COLUMNS,
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        for model_id in sorted(common_models):
            writer.writerow(model_metadata[model_id])

    crispr_count, crispr_missing = write_feature_matrix(
        CRISPR_FILE,
        OUTPUT_DIR / "crispr_effect_common_models.tsv",
        common_models,
        selected_gene_ids,
        crispr_columns,
        model_id_column=0,
    )
    expression_count, expression_missing = write_feature_matrix(
        EXPRESSION_FILE,
        OUTPUT_DIR / "expression_logp1_common_models.tsv",
        common_models,
        selected_gene_ids,
        expression_columns,
        model_id_column=3,
        expression=True,
    )

    gene_coverage = {
        "label_genes": len(label_genes),
        "genes_in_both_feature_matrices": len(selected_gene_ids),
        "genes_missing_from_either_matrix": len(label_genes) - len(selected_gene_ids),
    }

    summary_path = OUTPUT_DIR / "preprocessing_summary.txt"
    with summary_path.open("w", encoding="utf-8") as file:
        file.write("DepMap preprocessing summary\n")
        file.write(f"Input label file: {LABEL_FILE}\n")
        file.write(f"CRISPR source: {CRISPR_FILE}\n")
        file.write(f"Expression source: {EXPRESSION_FILE}\n")
        file.write(f"Model metadata source: {MODEL_FILE}\n")
        file.write(f"Common models: {len(common_models)}\n")
        file.write(f"Label genes: {gene_coverage['label_genes']}\n")
        file.write(
            "Genes present in both feature matrices: "
            f"{gene_coverage['genes_in_both_feature_matrices']}\n"
        )
        file.write(
            "Label genes missing from at least one feature matrix: "
            f"{gene_coverage['genes_missing_from_either_matrix']}\n"
        )
        file.write(
            f"CRISPR rows written: {crispr_count}; missing numeric values: "
            f"{crispr_missing}\n"
        )
        file.write(
            f"Expression rows written: {expression_count}; missing numeric "
            f"values: {expression_missing}\n"
        )
        file.write(
            "Expression filtering: IsDefaultEntryForModel == Yes\n"
        )
        file.write(
            "Feature scaling/imputation: not applied; fit these on training "
            "data only in a later modeling step.\n"
        )

    print(f"Model chung: {len(common_models)}")
    print(f"Gene trong bảng nhãn: {len(label_genes)}")
    print(f"Gene có trong cả hai ma trận: {len(selected_gene_ids)}")
    print(
        "Gene thiếu ở ít nhất một ma trận: "
        f"{gene_coverage['genes_missing_from_either_matrix']}"
    )
    print(
        f"CRISPR: {crispr_count} models, "
        f"{crispr_missing} giá trị thiếu/không hữu hạn"
    )
    print(
        f"Expression: {expression_count} models, "
        f"{expression_missing} giá trị thiếu/không hữu hạn"
    )
    print(f"Kết quả: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
