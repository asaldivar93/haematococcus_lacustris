
"""
Parse CDS features from a GenBank file into CSV.

Usage:
    python genbank_cds_to_csv.py input.gbk output.csv

Install dependency:
    pip install biopython
"""

import argparse
import csv
import json
from pathlib import Path
from Bio import SeqIO


def get_first(qualifiers, key, default=""):
    """Return the first value for a GenBank qualifier."""
    values = qualifiers.get(key, [])
    if not values:
        return default
    return values[0]


def parse_cds_features(genbank_file):
    """Yield dictionaries for each CDS feature in a GenBank file."""
    for record in SeqIO.parse(genbank_file, "genbank"):
        for feature in record.features:
            if feature.type != "CDS":
                continue

            qualifiers = feature.qualifiers

            locus_tag = get_first(qualifiers, "locus_tag")
            protein_id = get_first(qualifiers, "protein_id")
            gene = get_first(qualifiers, "gene")
            product = get_first(qualifiers, "product")
            note = get_first(qualifiers, "note")
            function = get_first(qualifiers, "function")
            ec_number = "; ".join(qualifiers.get("EC_number", []))
            db_xref = "; ".join(qualifiers.get("db_xref", []))
            translation = get_first(qualifiers, "translation")

            row = {
                "record_id": record.id,
                #"record_name": record.name,
                #"record_description": record.description,
                "locus_tag": locus_tag,
                "protein_id": protein_id,
                "gene": gene,
                "name": product,
                #"product": product,
                #"note": note,
                "function": function,
                "ec_number": ec_number,
                "db_xref": db_xref,
                #"start": int(feature.location.start) + 1,  # GenBank-style 1-based start
                #"end": int(feature.location.end),
                #"strand": feature.location.strand,
                #"location": str(feature.location),
                #"translation": translation,
                #"all_cds_annotations_json": json.dumps(qualifiers, ensure_ascii=False),
            }

            yield row


def write_csv(rows, output_csv):
    fieldnames = [
        "record_id",
        #"record_name",
        #"record_description",
        "locus_tag",
        "protein_id",
        "gene",
        "name",
        #"product",
        #"note",
        "function",
        "ec_number",
        "db_xref",
        #"start",
        #"end",
        #"strand",
        #"location",
        #"translation",
        "all_cds_annotations_json",
    ]

    with open(output_csv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()

        for row in rows:
            writer.writerow(row)


def main():
    parser = argparse.ArgumentParser(
        description="Extract CDS features from a GenBank file into CSV."
    )
    parser.add_argument("genbank_file", help="Input GenBank file, e.g. genome.gbk")
    parser.add_argument("output_csv", help="Output CSV file, e.g. cds_features.csv")

    args = parser.parse_args()

    genbank_path = Path(args.genbank_file)
    if not genbank_path.exists():
        raise FileNotFoundError(f"Input file not found: {genbank_path}")

    rows = parse_cds_features(genbank_path)
    write_csv(rows, args.output_csv)

    print(f"Finished writing CDS annotations to: {args.output_csv}")


if __name__ == "__main__":
    main()
