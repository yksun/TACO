"""Final report generation for TACO."""
import os
import sys
import argparse
import csv
from collections import defaultdict


def read_metric_merged(csv_path):
    """Read a merged metric CSV file.

    Args:
        csv_path: Path to merged metric CSV file

    Returns:
        dict: Metric data with appropriate structure
    """
    if not os.path.exists(csv_path):
        return {}

    data = {}
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row and any(row.values()):
                # Store row with all fields
                for key, value in row.items():
                    if key:
                        data[key] = value

    return data


def read_assembly_info(csv_path):
    """Read assembly_info.csv and return list of row dicts.

    Args:
        csv_path: Path to assembly_info.csv

    Returns:
        list: List of row dictionaries
    """
    rows = []

    if not os.path.exists(csv_path):
        return rows

    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row and any(row.values()):
                rows.append(row)

    return rows


def generate_final_report(output_path, assembly_info_path, decision_file=None):
    """Generate final combined report.

    Args:
        output_path: Output final_result.csv path
        assembly_info_path: Path to assembly_info.csv
        decision_file: Optional decision text file path
    """
    # Read assembly info
    assembly_rows = read_assembly_info(assembly_info_path)

    # Read decision if provided
    decision = None
    if decision_file and os.path.exists(decision_file):
        with open(decision_file, 'r') as f:
            decision = f.read().strip()

    # Build final report
    final_rows = []

    # Add assembly rows
    for row in assembly_rows:
        final_rows.append(row)

    # Add decision if available
    if decision:
        decision_row = {"metric": "selected_backbone", "value": decision}
        final_rows.append(decision_row)

    # Write final report
    if final_rows:
        fieldnames = set()
        for row in final_rows:
            fieldnames.update(row.keys())
        fieldnames = sorted(list(fieldnames))

        with open(output_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(final_rows)

        print(f"Final report written to {output_path}")
    else:
        print("No data to write in final report")


def main():
    """Command-line interface for report generation."""
    parser = argparse.ArgumentParser(description="TACO final report generation")
    parser.add_argument("--output", required=True, help="Output report path")
    parser.add_argument("--assembly-info", required=True, help="Path to assembly_info.csv")
    parser.add_argument("--decision-file", help="Optional decision text file")

    args = parser.parse_args()

    generate_final_report(
        args.output,
        args.assembly_info,
        decision_file=args.decision_file,
    )


if __name__ == "__main__":
    main()
