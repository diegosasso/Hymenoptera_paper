"""Extract published HAIKA matrices without changing source values or taxon names.

Run from the repository root with Python, openpyxl and numpy. This is lightweight
source-data extraction; evolutionary-model fitting belongs on Roihu.
"""

from pathlib import Path
import csv
import hashlib
import io
import json
import shutil
import zipfile

import numpy as np
import openpyxl


source = Path("/Users/taravser/Library/CloudStorage/OneDrive-UniversityofHelsinki/My_papers/Hymenoptera-Onto-Phylo/HAIKA-wing-chars/HAIKA-paper")
root = Path("comparisons/haika-wings")
raw = root / "data-raw"
data = root / "data"
provenance = root / "provenance"
raw.mkdir(parents=True, exist_ok=True)
data.mkdir(parents=True, exist_ok=True)
provenance.mkdir(parents=True, exist_ok=True)

source_files = {
    "haika-raw-matrices.xlsx": source / "S3. Raw matrices used in case studies 1–3 (1).xlsx",
    "haika-case2-phylogenetics.zip": source / "S4. Phylogenetic matrix and results for Hymenoptera in Case 2 (1).zip",
    "haika-supplement.pdf": source / "sciadv.aef0605_sm.pdf",
    "haika-supplement-tables.zip": source / "sciadv.aef0605_tables_s2_to_s11.zip",
}
checksums = []
for destination, original in source_files.items():
    target = raw / destination
    shutil.copy2(original, target)
    checksums.append([str(original), str(target), original.stat().st_size,
                      hashlib.sha256(original.read_bytes()).hexdigest()])

with (provenance / "haika-source-checksums.csv").open("w", newline="") as stream:
    writer = csv.writer(stream)
    writer.writerow(["source_file", "archived_file", "bytes", "sha256"])
    writer.writerows(checksums)

workbook = openpyxl.load_workbook(raw / "haika-raw-matrices.xlsx", read_only=True, data_only=True)
summary = {}
all_metadata = []
case_matrices = {}
for sheet_name, case_name in [("Case1", "case1"), ("Case 2", "case2"), ("Case3", "case3")]:
    worksheet = workbook[sheet_name]
    records = []
    metadata = []
    for row_number, row in enumerate(worksheet.values, start=1):
        if row_number < 3 or row[0] is None:
            continue
        label = str(row[0])
        group = str(row[1])
        values = list(row[2:])
        assert len(values) == 28 and all(isinstance(value, (int, float)) for value in values)
        records.append([label] + values)
        metadata.append([case_name, label, group, sheet_name, row_number])
    feature_names = ["feature_%02d" % index for index in range(1, 29)]
    with (data / ("haika-" + case_name + "-features.csv")).open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["taxon_label"] + feature_names)
        writer.writerows(records)
    all_metadata.extend(metadata)
    matrix = np.array([row[1:] for row in records], dtype=float)
    case_matrices[case_name] = (records, matrix)
    unique_vectors, duplicate_id, counts = np.unique(matrix, axis=0, return_inverse=True, return_counts=True)
    with (data / ("haika-" + case_name + "-duplicates.csv")).open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["taxon_label", "identical_vector_group", "number_with_identical_vector"])
        for row, identity in zip(records, duplicate_id):
            writer.writerow([row[0], int(identity) + 1, int(counts[identity])])
    summary[case_name] = {
        "source_sheet": sheet_name, "rows": len(records), "features": matrix.shape[1],
        "unique_taxon_labels": len(set(row[0] for row in records)),
        "unique_feature_vectors": len(unique_vectors),
        "constant_features": int(np.sum(np.ptp(matrix, axis=0) == 0)),
        "minimum": float(matrix.min()), "maximum": float(matrix.max()),
    }

with (data / "haika-taxonomy-source.csv").open("w", newline="") as stream:
    writer = csv.writer(stream)
    writer.writerow(["dataset", "taxon_label", "published_group", "source_sheet", "source_row"])
    writer.writerows(all_metadata)

with zipfile.ZipFile(raw / "haika-supplement-tables.zip") as archive:
    spreadsheet = archive.read("aef0605_tables_s2_to_s11.xlsx")
    (raw / "haika-supplement-tables.xlsx").write_bytes(spreadsheet)
tables = openpyxl.load_workbook(io.BytesIO(spreadsheet), read_only=True, data_only=True)
pc_records = []
for row_number, row in enumerate(tables["S5"].values, start=1):
    if row_number >= 4 and row[0] is not None:
        assert isinstance(row[1], (int, float)) and isinstance(row[2], (int, float))
        pc_records.append([str(row[0]), float(row[1]), float(row[2]), row_number])
with (data / "haika-case3-published-pc-scores.csv").open("w", newline="") as stream:
    writer = csv.writer(stream)
    writer.writerow(["taxon_label", "PC1", "PC2", "source_row_S5"])
    writer.writerows(pc_records)

# Preserve the published 64-dimensional matrix separately; it is the Case 2
# phylogenetic matrix and must not be concatenated with the 28-column Case 3 data.
with zipfile.ZipFile(raw / "haika-case2-phylogenetics.zip") as archive:
    member = "S4. Phylogenetic matrix and results for Hymenoptera in Case 2/matrix_Continuous (HAIKA).txt"
    continuous_text = archive.read(member).decode("utf-8-sig")
    (raw / "haika-case2-continuous64.txt").write_text(continuous_text)
continuous_records = []
for line in continuous_text.splitlines():
    if not line.strip():
        continue
    parts = line.split()
    assert len(parts) == 65
    continuous_records.append([parts[0]] + [float(value) for value in parts[1:]])
with (data / "haika-case2-phylogenetic64.csv").open("w", newline="") as stream:
    writer = csv.writer(stream)
    writer.writerow(["taxon_label"] + ["feature_%02d" % index for index in range(1, 65)])
    writer.writerows(continuous_records)
matrix64 = np.array([row[1:] for row in continuous_records])
summary["case2_phylogenetic64"] = {
    "rows": len(continuous_records), "features": 64,
    "unique_feature_vectors": len(np.unique(matrix64, axis=0)),
    "constant_features": int(np.sum(np.ptp(matrix64, axis=0) == 0)),
}

# Compare source-matched 28-feature PCA to published Case 3 scores. This is a
# provenance diagnostic, not a replacement for the publication's original axes.
case3_records, case3_matrix = case_matrices["case3"]
case3_labels = [row[0] for row in case3_records]
pc_lookup = {row[0]: row for row in pc_records}
assert len(pc_lookup) == len(pc_records)
assert set(case3_labels) == set(pc_lookup)
published_pc = np.array([pc_lookup[label][1:3] for label in case3_labels])
diagnostic_rows = []
for scaling in ["center_only", "center_and_unit_variance"]:
    centered = case3_matrix - case3_matrix.mean(axis=0)
    if scaling == "center_and_unit_variance":
        centered = centered / case3_matrix.std(axis=0, ddof=1)
    u, singular, vt = np.linalg.svd(centered, full_matrices=False)
    scores = u * singular
    explained = singular ** 2 / np.sum(singular ** 2)
    for axis in range(2):
        correlation = np.corrcoef(scores[:, axis], published_pc[:, axis])[0, 1]
        diagnostic_rows.append([scaling, axis + 1, explained[axis], correlation])
with (provenance / "haika-pca-reproduction-diagnostic.csv").open("w", newline="") as stream:
    writer = csv.writer(stream)
    writer.writerow(["feature_scaling", "axis", "variance_fraction", "correlation_with_published_axis_sign_arbitrary"])
    writer.writerows(diagnostic_rows)

summary["published_pc_scores"] = {
    "rows": len(pc_records), "unique_pc_pairs": len(np.unique(published_pc, axis=0)),
    "taxon_sets_identical_to_case3": True,
}
(provenance / "haika-data-summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2))

# Map identical published vectors back to the broader Case 2 image-label matrix.
# Equality demonstrates reuse of a rounded vector, not specimen identity.
case2_records, case2_matrix = case_matrices["case2"]
case2_by_vector = {}
for record in case2_records:
    case2_by_vector.setdefault(tuple(record[1:]), []).append(record[0])
with (data / "haika-case3-vector-provenance.csv").open("w", newline="") as stream:
    writer = csv.writer(stream)
    writer.writerow(["case3_species", "n_identical_case2_rows", "identical_case2_labels"])
    for record in case3_records:
        labels = case2_by_vector.get(tuple(record[1:]), [])
        writer.writerow([record[0], len(labels), ";".join(labels)])

# All five independent-seed matrices are archived 64D image embeddings, with
# tab-separated image identifiers (which can themselves contain spaces).
seed_summary = []
with zipfile.ZipFile(raw / "haika-five-seeds.zip") as archive:
    for member in archive.namelist():
        records = []
        for line in archive.read(member).decode("utf-8-sig").splitlines():
            label, values = line.split("\t")
            features = [float(value) for value in values.split()]
            assert len(features) == 64
            records.append([label] + features)
        seed = member.removeprefix("image_embeddings_seed").removesuffix(".txt")
        with (data / ("haika-case2-seed" + seed + "-features64.csv")).open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(["taxon_label"] + ["feature_%02d" % index for index in range(1, 65)])
            writer.writerows(records)
        seed_matrix = np.array([row[1:] for row in records])
        seed_summary.append({"seed": seed, "rows": len(records),
                             "unique_vectors": len(np.unique(seed_matrix, axis=0)),
                             "constant_features": int(np.sum(np.ptp(seed_matrix, axis=0) == 0)),
                             "same_taxon_labels_as_case2": set(row[0] for row in records) == set(row[0] for row in case2_records)})
(provenance / "haika-seed-summary.json").write_text(json.dumps(seed_summary, indent=2) + "\n")

# Read only the archived code files needed to document rounding and PCA behavior.
with zipfile.ZipFile(raw / "haika-code.zip") as archive:
    for member in ["main.py", "README.md", "preprocessing.py"]:
        (provenance / ("haika-archived-" + member)).write_bytes(archive.read(member))

# Check the public archive downloads and supplied S3/S4 copies against the
# publisher's Zenodo checksums. The supplement workbook is a Science attachment.
zenodo = json.loads((provenance / "haika-zenodo-record.json").read_text())
zenodo_keys = {
    "haika-raw-matrices.xlsx": "S3. Raw matrices used in case studies 1–3.xlsx",
    "haika-case2-phylogenetics.zip": "S4. Phylogenetic matrix and results for Hymenoptera in Case 2.zip",
    "haika-five-seeds.zip": "S8. Embeddings 5random seeds.zip",
    "haika-code.zip": "S1. The code, sample files, and user guide for the HAIKA model.zip",
}
zenodo_files = {item["key"]: item for item in zenodo["files"]}
with (provenance / "haika-zenodo-checksums.csv").open("w", newline="") as stream:
    writer = csv.writer(stream)
    writer.writerow(["archived_file", "zenodo_key", "url", "bytes", "md5", "published_checksum", "matches_published_checksum", "sha256"])
    for local_name, archive_key in zenodo_keys.items():
        contents = (raw / local_name).read_bytes()
        file_metadata = zenodo_files[archive_key]
        md5 = hashlib.md5(contents).hexdigest()
        assert "md5:" + md5 == file_metadata["checksum"]
        writer.writerow([str(raw / local_name), archive_key, file_metadata["links"]["self"],
                         len(contents), md5, file_metadata["checksum"], True,
                         hashlib.sha256(contents).hexdigest()])
