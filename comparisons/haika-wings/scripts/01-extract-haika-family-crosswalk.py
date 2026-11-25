"""Build an explicit family-proxy crosswalk, authorized as a sensitivity analysis.

Only family names explicitly present in source image labels are accepted. No
taxonomic family is inferred from a source genus, tribe, or subfamily here.
"""

from pathlib import Path
import csv
import json
import re

root = Path("comparisons/haika-wings")
data = root / "data"
provenance = root / "provenance"
with (data / "tree-tip-metadata.csv").open() as stream:
    tips = list(csv.DictReader(stream))
with (data / "haika-case2-features.csv").open() as stream:
    case2 = list(csv.DictReader(stream))
with (data / "haika-case2-seed3407-features64.csv").open() as stream:
    seed = list(csv.DictReader(stream))

# These are literal differences between the archived S3 and S8 labels.
# The whole suffixes otherwise match one-to-one; preserve both original labels.
label_lookup = {}
label_records = []
for row in seed:
    original = row["taxon_label"]
    harmonized = original.replace("Apocrita_nn_Aculeata-", "Apocrita-")
    harmonized = harmonized.replace("Aculeata-Sphecidae Ammophilinae",
                                    "Aculeata-Sphecidae Arnmophilinae")
    assert harmonized not in label_lookup
    label_lookup[harmonized] = original
    label_records.append([harmonized, original,
                          "prefix_and_spelling_harmonization" if original != harmonized else "identical"])
assert set(label_lookup) == set(row["taxon_label"] for row in case2)
with (data / "haika-case2-seed-label-crosswalk.csv").open("w", newline="") as stream:
    writer = csv.writer(stream)
    writer.writerow(["case2_S3_label", "seed_S8_label", "mapping_basis"])
    writer.writerows(label_records)

rows = []
long_rows = []
family_tokens = {}
for row in case2:
    label = row["taxon_label"]
    family_tokens[label] = set(re.split(r"[^A-Za-z]+", label))

for tip in tips:
    morphology_family = tip["Family_Sharkey"]
    molecular_family = tip["Family_Blaimer"]
    morphology_matches = []
    molecular_matches = []
    for source_row in case2:
        label = source_row["taxon_label"]
        if morphology_family in family_tokens[label]:
            morphology_matches.append(label)
            long_rows.append([tip["tip_id"], tip["analysis_label"], morphology_family,
                              label, label_lookup[label], "explicit_family_token",
                              "morphology_family", "family_proxy_not_species_observation"])
        if molecular_family in family_tokens[label]:
            molecular_matches.append(label)
    rows.append({
        "tip_id": tip["tip_id"],
        "analysis_label": tip["analysis_label"],
        "molecular_species": tip["molecular_species"],
        "Family_Sharkey": morphology_family,
        "Family_Blaimer": molecular_family,
        "family_names_differ": morphology_family != molecular_family,
        "source_taxonomic_level": "family",
        "selected_family_basis": "Family_Sharkey",
        "selected_family": morphology_family,
        "n_source_images": len(morphology_matches),
        "source_S3_labels": ";".join(morphology_matches),
        "source_S8_labels": ";".join(label_lookup[label] for label in morphology_matches),
        "n_source_images_molecular_family": len(molecular_matches),
        "source_S3_labels_molecular_family": ";".join(molecular_matches),
        "mapping_status": "explicit_family_proxy" if morphology_matches else "unmatched_no_explicit_family_token",
        "assumption": "Source images represent the named family; no measured species or genus identity is implied.",
    })

with (data / "haika-family-proxy-crosswalk.csv").open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
with (data / "haika-family-proxy-source-rows.csv").open("w", newline="") as stream:
    writer = csv.writer(stream)
    writer.writerow(["tip_id", "analysis_label", "selected_family", "source_S3_label",
                     "source_S8_label", "mapping_basis", "target_family_basis", "interpretation"])
    writer.writerows(long_rows)

matched = [row for row in rows if row["n_source_images"] > 0]
summary = {
    "total_study_tips": len(rows),
    "matched_tips_Family_Sharkey": len(matched),
    "distinct_matched_families": len(set(row["selected_family"] for row in matched)),
    "distinct_source_images": len(set(row[3] for row in long_rows)),
    "matched_tips_Family_Blaimer": sum(row["n_source_images_molecular_family"] > 0 for row in rows),
    "unmatched": [{"tip_id": row["tip_id"], "genus": row["analysis_label"],
                    "family": row["selected_family"]} for row in rows if row["n_source_images"] == 0],
    "mapping_rule": "Exact explicit family tokens. No inferred taxonomy, synonym remapping, or species-level imputation.",
    "repeated_profile_warning": "All study tips in a family inherit the same family aggregate; they are not independent measured species."
}
(provenance / "haika-family-proxy-summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2))

# Exact genus-token overlap is reported separately and never promoted from a
# family placeholder such as the molecular label "Megalyridae".
genus_rows = []
for tip in tips:
    morphology_genus = tip["Genus_Sharkey"]
    molecular_genus = tip["Genus_Blaimer"]
    morphology_matches = []
    molecular_matches = []
    for source_row in case2:
        label = source_row["taxon_label"]
        if morphology_genus in family_tokens[label] and not morphology_genus.endswith("idae"):
            morphology_matches.append(label)
        if molecular_genus in family_tokens[label] and not molecular_genus.endswith("idae"):
            molecular_matches.append(label)
    genus_rows.append([tip["tip_id"], morphology_genus, molecular_genus,
                       ";".join(morphology_matches), ";".join(molecular_matches)])
with (data / "haika-case2-genus-overlap.csv").open("w", newline="") as stream:
    writer = csv.writer(stream)
    writer.writerow(["tip_id", "morphology_genus", "molecular_genus",
                     "case2_matches_morphology_genus", "case2_matches_molecular_genus"])
    writer.writerows(genus_rows)
