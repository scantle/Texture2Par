"""Apply a tiny deterministic X offset to co-located wells in a comparison dataset."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def jitter(input_path: Path, output_path: Path, spacing: float) -> tuple[int, int]:
    with input_path.open(newline="", encoding="utf-8") as source:
        reader = csv.DictReader(source)
        if reader.fieldnames is None:
            raise ValueError("Dataset has no header")
        rows = list(reader)
        fieldnames = reader.fieldnames

    well_column = "ID" if "ID" in fieldnames else "WellIndex"
    wells_by_xy: dict[tuple[str, str], set[str]] = defaultdict(set)
    for row in rows:
        wells_by_xy[(row["X"], row["Y"])].add(row[well_column])

    offsets: dict[tuple[str, str, str], float] = {}
    affected_groups = 0
    affected_wells = 0
    for (x, y), wells in wells_by_xy.items():
        ordered = sorted(wells, key=lambda value: int(value))
        if len(ordered) < 2:
            continue
        affected_groups += 1
        affected_wells += len(ordered)
        center = 0.5 * (len(ordered) - 1)
        for index, well in enumerate(ordered):
            offsets[(x, y, well)] = (index - center) * spacing

    for row in rows:
        key = (row["X"], row["Y"], row[well_column])
        if key in offsets:
            row["X"] = f"{float(row['X']) + offsets[key]:.10f}"

    with output_path.open("w", newline="", encoding="utf-8") as target:
        writer = csv.DictWriter(target, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    return affected_groups, affected_wells


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--spacing", type=float, default=0.01)
    args = parser.parse_args()
    groups, wells = jitter(args.input, args.output, args.spacing)
    print(f"jittered {wells} wells in {groups} co-located XY groups")


if __name__ == "__main__":
    main()
