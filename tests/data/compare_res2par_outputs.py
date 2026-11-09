"""Compare staged legacy and libkriging Res2Par CSV outputs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


NAMES = [
    "LOGRHO",
    "FINE",
    "MIXED_FINE",
    "SAND",
    "MIXED_COARSE",
    "VERY_COARSE",
    "KH",
    "KV",
    "SS",
    "SY",
]


def metrics(expected: np.ndarray, actual: np.ndarray) -> dict[str, float | int]:
    valid = np.isfinite(expected) & np.isfinite(actual) & (expected > -999) & (actual > -999)
    difference = actual[valid] - expected[valid]
    if difference.size == 0:
        return {"count": 0}
    correlation = float(np.corrcoef(expected[valid], actual[valid])[0, 1])
    absolute = np.abs(difference)
    return {
        "count": int(difference.size),
        "mean_abs": float(np.mean(absolute)),
        "rmse": float(np.sqrt(np.mean(difference**2))),
        "p95_abs": float(np.percentile(absolute, 95)),
        "p99_abs": float(np.percentile(absolute, 99)),
        "p99_9_abs": float(np.percentile(absolute, 99.9)),
        "max_abs": float(np.max(absolute)),
        "count_abs_gt_0_01": int(np.count_nonzero(absolute > 0.01)),
        "correlation": correlation,
    }


def compare(baseline: Path, candidate: Path) -> dict[str, object]:
    result: dict[str, object] = {}
    for name in NAMES:
        expected = pd.read_csv(baseline / f"t2p_{name}.csv")
        actual = pd.read_csv(candidate / f"svihm_{name}.csv")
        columns = [column for column in expected if column.startswith(("Layer", "Var_Layer"))]
        result[name] = {
            column: metrics(expected[column].to_numpy(float), actual[column].to_numpy(float))
            for column in columns
        }
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("baseline", type=Path)
    parser.add_argument("candidate", type=Path)
    args = parser.parse_args()
    print(json.dumps(compare(args.baseline, args.candidate), indent=2))


if __name__ == "__main__":
    main()
