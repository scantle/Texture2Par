"""Convert a legacy Texture2Par bottom-depth-only dataset for comparison runs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def convert(input_path: Path, output_path: Path) -> tuple[int, int]:
    written = 0
    skipped_zero_thickness = 0
    previous_bottom: dict[str, float] = {}

    with input_path.open(newline="", encoding="utf-8") as source:
        reader = csv.DictReader(source)
        if reader.fieldnames is None or "Depth" not in reader.fieldnames:
            raise ValueError("Legacy dataset must contain a Depth column")
        well_column = "ID" if "ID" in reader.fieldnames else "WellIndex"
        output_fields = [
            field for field in reader.fieldnames if field != "Depth"
        ]
        depth_index = reader.fieldnames.index("Depth")
        output_fields[depth_index:depth_index] = ["TopDepth", "BotDepth"]

        with output_path.open("w", newline="", encoding="utf-8") as target:
            writer = csv.DictWriter(target, fieldnames=output_fields)
            writer.writeheader()
            for row in reader:
                well = row[well_column]
                bottom = float(row["Depth"])
                top = previous_bottom.get(well, 0.0)
                if bottom < top:
                    raise ValueError(
                        f"Bottom depth decreases for well {well}: {top} -> {bottom}"
                    )
                previous_bottom[well] = bottom
                if bottom == top:
                    skipped_zero_thickness += 1
                    continue
                converted = {
                    field: value for field, value in row.items() if field != "Depth"
                }
                converted["TopDepth"] = f"{top:.10g}"
                converted["BotDepth"] = f"{bottom:.10g}"
                writer.writerow(converted)
                written += 1

    return written, skipped_zero_thickness


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    written, skipped = convert(args.input, args.output)
    print(f"wrote {written} rows; skipped {skipped} zero-thickness rows")


if __name__ == "__main__":
    main()
