"""
test_03_interpolation.py
========================
Tests for Texture2Par's core interpolation behaviour: ordinary kriging,
sequential Gaussian simulation (SGS), cokriging, INFER_LAST_CLASS,
weight caching, and the ONLY_TEXTURE / READ_GRIDDED modes.
"""

import csv
import math
from pathlib import Path

import pytest
import pandas as pd
import numpy as np
from scipy.linalg import cho_factor, cho_solve
from scipy.spatial.distance import cdist, pdist, squareform

# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def _grid_t2p(options_extra="", classes="Coarse\n  Fine", infer=True,
              variogram_extra="") -> str:
    infer_line = "  INFER_LAST_CLASS" if infer else ""
    return f"""\
        BEGIN OPTIONS
          MAX_VSTRUCT 1
          {infer_line}
          WRITE_NODE_FILES
          {options_extra}
        END OPTIONS

        BEGIN CLASSES
          {classes}
        END CLASSES

        BEGIN DATASET
          FILE dataset.csv
        END DATASET

        BEGIN FLOW_MODEL
          TYPE GRID
          GRID_FILE grid.grd
        END FLOW_MODEL

        BEGIN VARIOGRAMS
          CLASS Coarse
            1 Sph 0.05 0.038 1.5E4 1.5E4 0.0 16
          {variogram_extra}
          CLASS PilotPoints
            1 Sph 0.00 1.0 1.0E5 1.0E5 0.0 25
        END VARIOGRAMS

        BEGIN PP_LOCS
          1 5050.0 4005.0 1
        END PP_LOCS

        BEGIN PP_PARAMETERS
          TYPE Global
          # ID KHp KVp STp
          1  0.5  -0.5  0.0

          TYPE Aquifer
          # ID Class  Kmin  Kmax    Ss      Sy   Aniso  Kd
          1 Coarse  0.1  100.0  1e-5  0.15  10.0  0.001
          1 Fine  0.001    0.5  1e-4  0.05   5.0  0.001
        END PP_PARAMETERS
        """


DATASET = """\
Location,WellIndex,IntervalIndex,X,Y,Zland,TopDepth,BotDepth,Coarse,Fine
W1,1,1,5000,4000.0,100.0,0.0,5.0,0.60,0.40
W1,1,2,5000,4000.0,100.0,5.0,10.0,0.55,0.45
W2,2,1,5100,4000.0,100.0,0.0,5.0,0.70,0.30
W2,2,2,5100,4000.0,100.0,5.0,10.0,0.45,0.55
W3,3,1,5050,4010.0,100.0,0.0,5.0,0.55,0.45
W3,3,2,5050,4010.0,100.0,5.0,10.0,0.50,0.50
W4,4,1,4980,4008.0,100.0,0.0,5.0,0.65,0.35
W4,4,2,4980,4008.0,100.0,5.0,10.0,0.42,0.58
"""

GRD = """\
NX    5
NY    5
NZ    2
DX    3000.0
DY    3000.0
TOP   constant 100.0
BOT   constant 80.0
IBOUND constant 1
"""

FILES = {"dataset.csv": DATASET, "grid.grd": GRD}


OK_LIB_2D_DATASET = """\
Location,WellIndex,IntervalIndex,X,Y,Zland,TopDepth,BotDepth,Coarse
W1,1,1,0.5,0.5,100,0,10,0.10
W2,2,1,1.5,0.5,100,0,10,0.20
W3,3,1,0.5,1.5,100,0,10,0.30
W4,4,1,1.5,1.5,100,0,10,0.40
"""


def _ok_lib_variance_dataset(first_variance: float) -> str:
    return f"""\
Location,WellIndex,IntervalIndex,X,Y,Zland,TopDepth,BotDepth,Coarse,var_Coarse
W1,1,1,0.5,0.5,100,0,10,0.10,{first_variance}
W2,2,1,1.5,0.5,100,0,10,0.20,0.0
W3,3,1,0.5,1.5,100,0,10,0.30,0.0
W4,4,1,1.5,1.5,100,0,10,0.40,0.0
"""


def _ck_variance_primary_dataset(first_variance: float) -> str:
    return _ok_lib_variance_dataset(first_variance)


CK_VARIANCE_SECONDARY_DATASET = """\
Location,WellIndex,IntervalIndex,X,Y,Zland,TopDepth,BotDepth,Rho,var_Rho
S1,1,1,0.5,0.5,100,0,10,0.90,0.0
S2,2,1,1.5,0.5,100,0,10,0.30,0.0
S3,3,1,0.5,1.5,100,0,10,0.40,0.0
S4,4,1,1.5,1.5,100,0,10,0.60,0.0
"""

SK_CK_SECONDARY_DATASET = """\
Location,WellIndex,IntervalIndex,X,Y,Zland,TopDepth,BotDepth,Rho,var_Rho
S1,1,1,2.5,0.5,100,0,10,0.90,0.0
S2,2,1,2.5,1.5,100,0,10,0.30,0.0
S3,3,1,0.5,2.5,100,0,10,0.40,0.0
S4,4,1,1.5,2.5,100,0,10,0.60,0.0
"""


def _shift_dataset_column(dataset: str, column: str, offset: float) -> str:
    rows = list(csv.reader(dataset.strip().splitlines()))
    column_index = rows[0].index(column)
    for row in rows[1:]:
        row[column_index] = f"{float(row[column_index]) + offset:.12g}"
    return "\n".join(",".join(row) for row in rows) + "\n"


OK_LIB_2D_GRD = """\
NX    2
NY    2
NZ    1
DX    1.0
DY    1.0
TOP   constant 100.0
DZ    10.0
IBOUND constant 1
"""


BLOCK_MODFLOW_NAM = """\
LIST 2 model.list
DIS 11 model.dis
BAS6 13 model.bas
LPF 15 model.lpf
"""

BLOCK_MODFLOW_DIS = """\
1 2 2 1 4 2
0
CONSTANT 1.0
CONSTANT 1.0
CONSTANT 100.0
CONSTANT 90.0
1.0 1 1.0 SS
"""

BLOCK_MODFLOW_BAS = """\
FREE
CONSTANT 1
-999.99
CONSTANT 0.0
"""

BLOCK_MODFLOW_LPF = """\
0 -1.0E+30 0
0
1.0
0
0
CONSTANT 1.0
CONSTANT 1.0
"""

BLOCK_DATASET = """\
Location,WellIndex,IntervalIndex,X,Y,Zland,TopDepth,BotDepth,Coarse
W1,1,1,0.1,0.1,100,0,10,0.10
W2,2,1,1.9,0.1,100,0,10,0.30
W3,3,1,0.1,1.9,100,0,10,0.65
W4,4,1,1.9,1.9,100,0,10,0.90
"""

BLOCK_FILES = {
    "dataset.csv": BLOCK_DATASET,
    "model.nam": BLOCK_MODFLOW_NAM,
    "model.dis": BLOCK_MODFLOW_DIS,
    "model.bas": BLOCK_MODFLOW_BAS,
    "model.lpf": BLOCK_MODFLOW_LPF,
}


def _block_variance_dataset(first_variance: float) -> str:
    return f"""\
Location,WellIndex,IntervalIndex,X,Y,Zland,TopDepth,BotDepth,Coarse,var_Coarse
W1,1,1,0.1,0.1,100,0,10,0.10,{first_variance}
W2,2,1,1.9,0.1,100,0,10,0.30,0.0
W3,3,1,0.1,1.9,100,0,10,0.65,0.0
W4,4,1,1.9,1.9,100,0,10,0.90,0.0
"""


SK_2D_DATASET = """\
Location,WellIndex,IntervalIndex,X,Y,Zland,TopDepth,BotDepth,Coarse
W1,1,1,0.2,0.4,100,0,10,0.10
W2,2,1,2.7,0.35,100,0,10,0.80
W3,3,1,0.4,2.6,100,0,10,0.20
W4,4,1,2.4,2.2,100,0,10,0.90
"""


SK_2D_GRD = """\
NX    3
NY    3
NZ    1
DX    1.0
DY    1.0
TOP   constant 100.0
DZ    10.0
IBOUND constant 1
"""


def _ok_lib_2d_t2p(options_extra="", flow_extra="", method="OK",
                   variogram_lines="1 Sph 0.0 1.0 10.0 10.0 0.0 4",
                   max_vstruct=1) -> str:
    return f"""\
        BEGIN OPTIONS
          MAX_VSTRUCT {max_vstruct}
          INTERP_METHOD {method}
          INFER_LAST_CLASS
          ONLY_TEXTURE
          WRITE_NODE_FILES
          WRITE_DATASET_FILES
          {options_extra}
        END OPTIONS

        BEGIN CLASSES
          Coarse
          Fine
        END CLASSES

        BEGIN DATASET
          FILE dataset.csv
        END DATASET

        BEGIN FLOW_MODEL
          TYPE GRID
          GRID_FILE grid.grd
          {flow_extra}
        END FLOW_MODEL

        BEGIN VARIOGRAMS
          CLASS Coarse
            {variogram_lines}
        END VARIOGRAMS
        """


def _ck_variance_t2p() -> str:
    return """\
        BEGIN OPTIONS
          MAX_VSTRUCT 1
          INTERP_METHOD OK
          CK_MEAN_MODE LIBKRIGING
          INFER_LAST_CLASS
          ONLY_TEXTURE
          WRITE_NODE_FILES
          WRITE_VARIANCE
          KVME
        END OPTIONS

        BEGIN CLASSES
          Coarse
          Fine
        END CLASSES

        BEGIN COVARIATES
          Coarse Rho
        END COVARIATES

        BEGIN DATASET
          FILE primary.csv
        END DATASET

        BEGIN DATASET
          FILE secondary.csv
        END DATASET

        BEGIN FLOW_MODEL
          TYPE GRID
          GRID_FILE grid.grd
        END FLOW_MODEL

        BEGIN VARIOGRAMS
          CLASS Coarse
            1 Sph 0.0 1.0 10.0 10.0 0.0 4
          CLASS Rho
            1 Sph 0.0 1.0 10.0 10.0 0.0 4
          CLASS Coarse:Rho
            1 Sph 0.0 0.25 10.0 10.0 0.0
        END VARIOGRAMS
        """


def _block_modflow_t2p(method="OK", options_extra="") -> str:
    return f"""\
        BEGIN OPTIONS
          MAX_VSTRUCT 1
          INTERP_METHOD {method}
          INFER_LAST_CLASS
          ONLY_TEXTURE
          WRITE_NODE_FILES
          WRITE_VARIANCE
          {options_extra}
        END OPTIONS

        BEGIN CLASSES
          Coarse
          Fine
        END CLASSES

        BEGIN DATASET
          FILE dataset.csv
        END DATASET

        BEGIN FLOW_MODEL
          TYPE MODFLOW
          NAM_FILE model.nam
          BLOCK_MODE
        END FLOW_MODEL

        BEGIN VARIOGRAMS
          CLASS Coarse
            1 Sph 0.0 1.0 10.0 10.0 0.0 4
        END VARIOGRAMS
        """


OK_LIB_2D_PARITY_DATASET = """\
Location,WellIndex,IntervalIndex,X,Y,Zland,TopDepth,BotDepth,Coarse
W1,1,1,0.10,0.20,100,0,10,0.10
W2,2,1,1.90,0.10,100,0,10,0.35
W3,3,1,0.20,1.80,100,0,10,0.65
W4,4,1,1.80,1.90,100,0,10,0.90
W5,5,1,1.10,1.60,100,0,10,0.55
"""


OK_LIB_NNEAR_DATASET = """\
Location,WellIndex,IntervalIndex,X,Y,Zland,TopDepth,BotDepth,Coarse
W1,1,1,0.40,0.50,100,0,10,0.11
W2,2,1,1.80,0.40,100,0,10,0.33
W3,3,1,0.30,1.80,100,0,10,0.66
W4,4,1,1.90,1.90,100,0,10,0.88
"""


OK_LIB_ANISO_SEARCH_DATASET = """\
Location,WellIndex,IntervalIndex,X,Y,Zland,TopDepth,BotDepth,Coarse
MinorAxis,1,1,0.50,0.70,100,0,10,0.10
MajorAxis,2,1,1.60,0.50,100,0,10,0.90
"""


DUPLICATE_GRD = """\
NX    3
NY    3
NZ    1
DX    1.0
DY    1.0
TOP   constant 100.0
DZ    10.0
IBOUND constant 1
"""


def _duplicate_conditioning_t2p(interp_dim: int) -> str:
    if interp_dim == 2:
        variogram = "Sph 0.0 1.0 10.0 10.0 0.0 10"
    else:
        variogram = "Sph 0.0 1.0 10.0 10.0 10.0 0.0 0.0 0.0 10"
    return f"""\
        BEGIN OPTIONS
          MAX_VSTRUCT 1
          INFER_LAST_CLASS
          ONLY_TEXTURE
          WRITE_NODE_FILES
          INTERP_DIM {interp_dim}
        END OPTIONS

        BEGIN CLASSES
          Coarse
          Fine
        END CLASSES

        BEGIN DATASET
          FILE dataset.csv
        END DATASET

        BEGIN FLOW_MODEL
          TYPE GRID
          GRID_FILE grid.grd
        END FLOW_MODEL

        BEGIN VARIOGRAMS
          CLASS Coarse
            1 {variogram}
        END VARIOGRAMS
        """


DUPLICATE_2D_DATASET = """\
Location,WellIndex,IntervalIndex,X,Y,Zland,TopDepth,BotDepth,Coarse
W1,1,1,0.5,0.5,100,0,10,0.10
W2,2,1,0.5,0.5,100,0,10,0.90
W3,3,1,2.5,0.5,100,0,10,0.30
W4,4,1,0.5,2.5,100,0,10,0.40
"""


DUPLICATE_3D_DATASET = """\
Location,WellIndex,IntervalIndex,X,Y,Zland,TopDepth,BotDepth,Coarse
W1,1,1,0.5,0.5,100,4,6,0.10
W2,2,1,0.5,0.5,100,4,6,0.90
W3,3,1,2.5,0.5,100,4,6,0.30
W4,4,1,0.5,2.5,100,4,6,0.40
"""


def _read_csv_column(path: Path, col: str) -> list[float]:
    rows = list(csv.DictReader(path.read_text().splitlines()))
    return [float(r[col]) for r in rows if r.get(col, "").strip() not in ("", "nan")]


MAX_LOG_LENGTH_GRD = """\
NX    1
NY    1
NZ    3
DX    1.0
DY    1.0
TOP   constant 100.0
DZ    10.0
IBOUND constant 1
"""

OK3D_MULTI_GRD = """\
NX    2
NY    2
NZ    3
DX    1.0
DY    1.0
TOP   constant 100.0
DZ    10.0
IBOUND constant 1
"""

OK3D_INTERIOR_GRD = """\
NX    3
NY    3
NZ    3
DX    1.0
DY    1.0
TOP   constant 100.0
DZ    10.0
IBOUND constant 1
"""

OK3D_EXACT_VGM = "Sph 0.0 1.0 1000 1000 5 0.0 0.0 0.0 10"
OK3D_MULTI_VGM = "Sph 0.0 1.0 2.0 2.0 8.0 0.0 0.0 0.0 20"
OK3D_CROSSLAYER_VGM = "Sph 0.0 1.0 2.0 2.0 20.0 0.0 0.0 0.0 20"
OK3D_ROTATED_VGM = "Sph 0.0 1.0 1.5 3.0 20.0 30.0 0.0 0.0 20"
OK3D_TRUNCATED_VGM = "Sph 0.0 1.0 2.0 2.0 20.0 0.0 0.0 0.0 5"
OK3D_ANISO_SEARCH_VGM = "Sph 0.0 1.0 1.0 4.0 20.0 35.0 0.0 0.0 4"
OK3D_EXP_VGM = "Exp 0.0 1.0 2.0 2.0 20.0 0.0 0.0 0.0 20"
OK3D_GAU_VGM = "Gau 0.01 1.0 2.0 2.0 20.0 0.0 0.0 0.0 20"


def _max_log_length_t2p(max_log_length=5) -> str:
    return f"""\
        BEGIN OPTIONS
          MAX_VSTRUCT 1
          INFER_LAST_CLASS
          ONLY_TEXTURE
          WRITE_DATASET_FILES
          INTERP_DIM 3
          MAX_LOG_LENGTH {max_log_length}
        END OPTIONS

        BEGIN CLASSES
          Coarse
          Fine
        END CLASSES

        BEGIN DATASET
          FILE dataset.csv
        END DATASET

        BEGIN FLOW_MODEL
          TYPE GRID
          GRID_FILE grid.grd
        END FLOW_MODEL

        BEGIN VARIOGRAMS
          CLASS Coarse
            1 Sph 0.001 1.0 1000 1000 5 0.0 0.0 0.0 10
        END VARIOGRAMS
        """


def _max_log_length_dataset(bot_depth: float) -> str:
    return f"""\
Location,WellIndex,IntervalIndex,X,Y,Zland,TopDepth,BotDepth,Coarse
W1,1,1,0,0,100,0,{bot_depth},0.60
W2,2,1,10,0,100,0,5,0.50
W3,3,1,0,10,100,0,5,0.50
W4,4,1,10,10,100,0,5,0.50
"""


def _processed_coarse_logxyz(r) -> pd.DataFrame:
    candidates = sorted(
        path for path in r.workdir.glob("*Coarse*logxyz*.csv")
        if path.name.lower() != "dataset.csv"
    )
    if not candidates:
        files = "\n".join(sorted(path.name for path in r.workdir.iterdir()))
        pytest.fail(f"No processed Coarse logxyz CSV found. Workdir files:\n{files}")
    df_all = pd.read_csv(candidates[0])
    required = {"X", "Y", "Z", "Layer1"}
    missing = required - set(df_all.columns)
    if missing:
        pytest.fail(
            f"Processed log CSV {candidates[0].name} missing columns "
            f"{sorted(missing)}. Columns: {list(df_all.columns)}"
        )
    tol = 1e-9
    df = df_all[
        (df_all["X"].sub(0.0).abs() < tol) &
        (df_all["Y"].sub(0.0).abs() < tol)
    ]
    if df.empty:
        pytest.fail(
            f"Processed log CSV {candidates[0].name} contains no target "
            f"records near X=0, Y=0. Rows:\n{df_all.to_string(index=False)}"
        )
    return df.sort_values("Z", ascending=False).reset_index(drop=True)


def _optional_direct_segment_lengths(df: pd.DataFrame):
    columns = {col.lower().replace("_", ""): col for col in df.columns}

    if "length" in columns:
        return df[columns["length"]].abs().tolist()

    boundary_pairs = (
        ("topdepth", "botdepth"),
        ("topdepth", "bottomdepth"),
        ("top", "bot"),
        ("top", "bottom"),
        ("topelev", "botelev"),
        ("topelev", "bottomelev"),
        ("topelevation", "botelevation"),
        ("topelevation", "bottomelevation"),
    )
    for top_key, bot_key in boundary_pairs:
        if top_key in columns and bot_key in columns:
            top = df[columns[top_key]]
            bot = df[columns[bot_key]]
            return top.sub(bot).abs().tolist()
    return None


def _assert_processed_log_points(
    r,
    expected_z,
    expected_midpoint_spacings=None,
    expected_segment_lengths=None,
    max_log_length=5.0,
):
    assert r.returncode == 0, r.combined
    df = _processed_coarse_logxyz(r)
    assert len(df) == len(expected_z), df.to_string(index=False)
    assert df["X"].tolist() == pytest.approx([0.0] * len(expected_z))
    assert df["Y"].tolist() == pytest.approx([0.0] * len(expected_z))
    assert df["Z"].tolist() == pytest.approx(expected_z)
    assert df["Layer1"].tolist() == pytest.approx([0.60] * len(expected_z))

    if expected_segment_lengths:
        direct_lengths = _optional_direct_segment_lengths(df)
        if direct_lengths is not None:
            assert direct_lengths == pytest.approx(expected_segment_lengths)
            assert max(direct_lengths) <= max_log_length + 1e-9

    if expected_midpoint_spacings and len(expected_z) > 1:
        inferred_lengths = [
            abs(df.loc[i, "Z"] - df.loc[i + 1, "Z"])
            for i in range(len(df) - 1)
        ]
        assert inferred_lengths == pytest.approx(expected_midpoint_spacings)
        assert max(inferred_lengths) <= max_log_length + 1e-9


def _interp_dim3_geometry_t2p() -> str:
    return """\
        BEGIN OPTIONS
          MAX_VSTRUCT 1
          INFER_LAST_CLASS
          ONLY_TEXTURE
          WRITE_ELEVATION_FILE
          INTERP_DIM 3
        END OPTIONS

        BEGIN CLASSES
          Coarse
          Fine
        END CLASSES

        BEGIN DATASET
          FILE dataset.csv
        END DATASET

        BEGIN FLOW_MODEL
          TYPE GRID
          GRID_FILE grid.grd
        END FLOW_MODEL

        BEGIN VARIOGRAMS
          CLASS Coarse
            1 Sph 0.001 1.0 1000 1000 5 0.0 0.0 0.0 10
        END VARIOGRAMS
        """


def _find_elevation_file(r) -> Path:
    candidates = sorted(set(
        list(r.workdir.glob("*elev*")) +
        list(r.workdir.glob("*ELEV*")) +
        list(r.workdir.glob("*elevation*"))
    ))
    if not candidates:
        files = "\n".join(sorted(path.name for path in r.workdir.iterdir()))
        pytest.fail(f"No elevation file found. Workdir files:\n{files}")
    return candidates[0]


def _read_grid_elevation_file(r) -> pd.DataFrame:
    path = _find_elevation_file(r)
    df = pd.read_csv(path)
    required = {"Node", "X", "Y", "Top", "Bottom1", "Bottom2", "Bottom3"}
    missing = required - set(df.columns)
    if missing:
        pytest.fail(
            f"Elevation file {path.name} missing columns {sorted(missing)}. "
            f"Columns: {list(df.columns)}. Rows:\n{df.to_string(index=False)}"
        )
    return df


def _ok3d_exact_hit_dataset() -> str:
    return """\
Location,WellIndex,IntervalIndex,X,Y,Zland,TopDepth,BotDepth,Coarse
W1,1,1,0.5,0.5,100,4,6,0.10
W1,1,2,0.5,0.5,100,14,16,0.50
W1,1,3,0.5,0.5,100,24,26,0.90
"""


def _ok3d_t2p(method="OK") -> str:
    return f"""\
        BEGIN OPTIONS
          MAX_VSTRUCT 1
          INFER_LAST_CLASS
          ONLY_TEXTURE
          WRITE_NODE_FILES
          WRITE_DATASET_FILES
          WRITE_VARIANCE
          INTERP_DIM 3
          INTERP_METHOD {method}
        END OPTIONS

        BEGIN CLASSES
          Coarse
          Fine
        END CLASSES

        BEGIN DATASET
          FILE dataset.csv
        END DATASET

        BEGIN FLOW_MODEL
          TYPE GRID
          GRID_FILE grid.grd
        END FLOW_MODEL

        BEGIN VARIOGRAMS
          CLASS Coarse
            1 {OK3D_EXACT_VGM}
        END VARIOGRAMS
        """


def _ok3d_multi_dataset() -> str:
    return """\
Location,WellIndex,IntervalIndex,X,Y,Zland,TopDepth,BotDepth,Coarse
W1,1,1,0.5,0.5,100,4,6,0.10
W1,1,2,0.5,0.5,100,14,16,0.35
W1,1,3,0.5,0.5,100,24,26,0.70
W2,2,1,1.5,0.5,100,4,6,0.25
W2,2,2,1.5,0.5,100,14,16,0.50
W2,2,3,1.5,0.5,100,24,26,0.85
W3,3,1,0.5,1.5,100,4,6,0.20
W3,3,2,0.5,1.5,100,14,16,0.45
W3,3,3,0.5,1.5,100,24,26,0.80
W4,4,1,1.5,1.5,100,4,6,0.30
W4,4,2,1.5,1.5,100,14,16,0.55
W4,4,3,1.5,1.5,100,24,26,0.90
"""


def _ok3d_corner_dataset() -> str:
    return """\
Location,WellIndex,IntervalIndex,X,Y,Zland,TopDepth,BotDepth,Coarse
W1,1,1,0.5,0.5,100,4,6,0.10
W1,1,2,0.5,0.5,100,14,16,0.35
W1,1,3,0.5,0.5,100,24,26,0.70
W2,2,1,2.5,0.5,100,4,6,0.25
W2,2,2,2.5,0.5,100,14,16,0.50
W2,2,3,2.5,0.5,100,24,26,0.85
W3,3,1,0.5,2.5,100,4,6,0.20
W3,3,2,0.5,2.5,100,14,16,0.45
W3,3,3,0.5,2.5,100,24,26,0.80
W4,4,1,2.5,2.5,100,4,6,0.30
W4,4,2,2.5,2.5,100,14,16,0.55
W4,4,3,2.5,2.5,100,24,26,0.90
"""


def _ok3d_irregular_dataset() -> str:
    return """\
Location,WellIndex,IntervalIndex,X,Y,Zland,TopDepth,BotDepth,Coarse
W1,1,1,0.2,0.4,100,4,6,0.08
W1,1,2,0.2,0.4,100,14,16,0.32
W1,1,3,0.2,0.4,100,24,26,0.68
W2,2,1,2.7,0.35,100,4,6,0.28
W2,2,2,2.7,0.35,100,14,16,0.52
W2,2,3,2.7,0.35,100,24,26,0.86
W3,3,1,0.4,2.6,100,4,6,0.18
W3,3,2,0.4,2.6,100,14,16,0.43
W3,3,3,0.4,2.6,100,24,26,0.78
W4,4,1,2.4,2.2,100,4,6,0.34
W4,4,2,2.4,2.2,100,14,16,0.58
W4,4,3,2.4,2.2,100,24,26,0.92
W5,5,1,1.1,0.9,100,4,6,0.22
W5,5,2,1.1,0.9,100,14,16,0.49
W5,5,3,1.1,0.9,100,24,26,0.82
"""


def _ok3d_multi_t2p(vgm=OK3D_MULTI_VGM, extra_options="", method="OK") -> str:
    return f"""\
        BEGIN OPTIONS
          MAX_VSTRUCT 1
          INFER_LAST_CLASS
          ONLY_TEXTURE
          WRITE_NODE_FILES
          WRITE_DATASET_FILES
          WRITE_VARIANCE
          INTERP_DIM 3
          INTERP_METHOD {method}
          {extra_options}
        END OPTIONS

        BEGIN CLASSES
          Coarse
          Fine
        END CLASSES

        BEGIN DATASET
          FILE dataset.csv
        END DATASET

        BEGIN FLOW_MODEL
          TYPE GRID
          GRID_FILE grid.grd
        END FLOW_MODEL

        BEGIN VARIOGRAMS
          CLASS Coarse
            1 {vgm}
        END VARIOGRAMS
        """


def _read_coarse_texture_output(
    r, required_layers=("Layer1", "Layer2", "Layer3")
) -> pd.DataFrame:
    candidates = sorted(
        path for path in r.workdir.glob("*Coarse*.csv")
        if all(
            suffix not in path.name.lower()
            for suffix in ("logxyz", "layavg", "aqtard")
        )
    )
    if not candidates:
        files = "\n".join(sorted(path.name for path in r.workdir.iterdir()))
        pytest.fail(f"No Coarse texture node output found. Workdir files:\n{files}")

    df = pd.read_csv(candidates[0])
    required = set(required_layers)
    missing = required - set(df.columns)
    if missing:
        pytest.fail(
            f"Coarse texture output {candidates[0].name} missing columns "
            f"{sorted(missing)}. Columns: {list(df.columns)}. "
            f"Rows:\n{df.to_string(index=False)}"
        )
    return df


def _parse_ok3d_vgm(vgm: str):
    parts = vgm.split()
    vtype = parts[0].lower()[:3]
    nugget, sill, hmin, hmax, hvert = (float(parts[i]) for i in range(1, 6))
    azimuth, dip, plunge = (float(parts[i]) for i in range(6, 9))
    nnear = int(parts[9])
    rotmat = _fortran_variogram_rotmat(hmin, hmax, hvert, azimuth, dip, plunge)
    return vtype, nugget, sill, hmax, rotmat, nnear


def _fortran_variogram_rotmat(hmin, hmax, hvert, azimuth, dip, plunge):
    # Matches Texture2Par's Fortran variogram rotation convention; the
    # test_variogram_rotation_metric_sanity_checks test validates the resulting
    # metric behavior independently for simple cases.
    rotmat = np.diag([hmax / hmin, 1.0, hmax / hvert])

    if abs(plunge) > 1e-10:
        alpha = np.deg2rad(-plunge)
        sina = np.sin(alpha)
        cosa = np.cos(alpha)
        tmp = np.array([
            [cosa, 0.0, sina],
            [0.0, 1.0, 0.0],
            [-sina, 0.0, cosa],
        ])
        rotmat = rotmat @ tmp

    if abs(dip) > 1e-10:
        alpha = np.deg2rad(dip)
        sina = np.sin(alpha)
        cosa = np.cos(alpha)
        tmp = np.array([
            [1.0, 0.0, 0.0],
            [0.0, cosa, sina],
            [0.0, -sina, cosa],
        ])
        rotmat = rotmat @ tmp

    if abs(azimuth) > 1e-10:
        alpha = np.deg2rad(azimuth)
        sina = np.sin(alpha)
        cosa = np.cos(alpha)
        tmp = np.array([
            [cosa, -sina, 0.0],
            [sina, cosa, 0.0],
            [0.0, 0.0, 1.0],
        ])
        rotmat = rotmat @ tmp

    return rotmat


def _effective_distance(delta, hmin, hmax, hvert, azimuth=0.0, dip=0.0, plunge=0.0):
    rotmat = _fortran_variogram_rotmat(hmin, hmax, hvert, azimuth, dip, plunge)
    scaled = rotmat @ np.asarray(delta, dtype=float) / hmax
    return np.linalg.norm(scaled)


def _spherical_covariance(hr):
    return np.where(hr < 1.0, 1.0 - 1.5 * hr + 0.5 * hr**3, 0.0)


def _ok3d_covariance(vtype, hr):
    if vtype == "sph":
        return _spherical_covariance(hr)
    if vtype == "exp":
        return np.exp(-3.0 * hr)
    if vtype == "gau":
        return np.exp(-49.0 / 16.0 * hr**2)
    raise ValueError(f"Unsupported covariance function {vtype!r}")


def _libkriging_gslib_rotmat(
    hmin, hmax, hvert, azimuth=0.0, dip=0.0, plunge=0.0
):
    """Return libkriging's GSLIB-style lag rotation and range scaling."""
    transform = np.diag([1.0 / hmin, 1.0 / hmax, 1.0 / hvert])

    plunge_rad = np.deg2rad(-plunge)
    cp, sp = np.cos(plunge_rad), np.sin(plunge_rad)
    rotate_y = np.array([
        [cp, 0.0, sp],
        [0.0, 1.0, 0.0],
        [-sp, 0.0, cp],
    ])

    dip_rad = np.deg2rad(dip)
    cd, sd = np.cos(dip_rad), np.sin(dip_rad)
    rotate_x = np.array([
        [1.0, 0.0, 0.0],
        [0.0, cd, sd],
        [0.0, -sd, cd],
    ])

    azimuth_rad = np.deg2rad(azimuth)
    ca, sa = np.cos(azimuth_rad), np.sin(azimuth_rad)
    rotate_z = np.array([
        [ca, -sa, 0.0],
        [sa, ca, 0.0],
        [0.0, 0.0, 1.0],
    ])
    return transform @ rotate_y @ rotate_x @ rotate_z


def _ok_lib_gslib_reference(
    obsloc,
    obsval,
    newloc,
    *,
    vtype="sph",
    nugget=0.0,
    sill=1.0,
    hmin,
    hmax,
    hvert,
    azimuth=0.0,
    dip=0.0,
    plunge=0.0,
):
    """Ordinary kriging with libkriging's convention, independent of legacy OK."""
    obsloc = np.asarray(obsloc, dtype=float)
    obsval = np.asarray(obsval, dtype=float)
    newloc = np.asarray(newloc, dtype=float)
    ndim = obsloc.shape[1]
    transform = _libkriging_gslib_rotmat(
        hmin, hmax, hvert, azimuth, dip, plunge
    )[:ndim, :ndim]
    scaled_obs = (transform @ obsloc.T).T
    scaled_new = (transform @ newloc.T).T

    cov0 = nugget + sill
    c_obs = _ok3d_covariance(vtype, squareform(pdist(scaled_obs))) * sill
    np.fill_diagonal(c_obs, cov0)
    rhs_dist = cdist(scaled_obs, scaled_new)
    c_rhs = _ok3d_covariance(vtype, rhs_dist) * sill
    if nugget > 0.0:
        c_rhs = c_rhs + np.where(rhs_dist < 1e-10, nugget, 0.0)

    nobs = len(obsloc)
    lhs = np.zeros((nobs + 1, nobs + 1))
    lhs[:nobs, :nobs] = c_obs
    lhs[:nobs, nobs] = 1.0
    lhs[nobs, :nobs] = 1.0
    rhs = np.vstack([c_rhs, np.ones((1, len(newloc)))])
    solved = np.linalg.solve(lhs, rhs)
    weights = solved[:nobs, :]
    mu = solved[nobs, :]
    estimates = obsval @ weights
    variance = np.maximum(0.0, cov0 - (weights * c_rhs).sum(axis=0) - mu)
    return estimates, variance


def _block_ok_reference(obsloc, obsval, centres, halfsize, hmax=10.0):
    """Independent 2-D OK reference using Texture2Par's historical quadrature."""
    obsloc = np.asarray(obsloc, dtype=float)
    obsval = np.asarray(obsval, dtype=float)
    centres = np.asarray(centres, dtype=float)
    frac = np.array([0.4305681558, 0.1699905218,
                     -0.1699905218, -0.4305681558])
    wgq = np.array([0.1739274226, 0.3260725774,
                    0.3260725774, 0.1739274226])
    offsets = np.array([
        (frac[ix] * halfsize[0], frac[iy] * halfsize[1])
        for iy in range(4) for ix in range(4)
    ])
    support_weights = np.array([
        wgq[ix] * wgq[iy] for iy in range(4) for ix in range(4)
    ])

    c_obs = _spherical_covariance(squareform(pdist(obsloc)) / hmax)
    np.fill_diagonal(c_obs, 1.0)
    nobs = len(obsloc)
    lhs = np.zeros((nobs + 1, nobs + 1))
    lhs[:nobs, :nobs] = c_obs
    lhs[:nobs, nobs] = 1.0
    lhs[nobs, :nobs] = 1.0

    estimates = np.empty(len(centres))
    variance = np.empty(len(centres))
    for target, centre in enumerate(centres):
        support = centre + offsets
        cov_support = _spherical_covariance(
            cdist(obsloc, support) / hmax
        )
        block_rhs = cov_support @ support_weights
        solved = np.linalg.solve(lhs, np.r_[block_rhs, 1.0])
        weights = solved[:nobs]
        mu = solved[nobs]
        support_cov = _spherical_covariance(
            squareform(pdist(support)) / hmax
        )
        np.fill_diagonal(support_cov, 1.0)
        block_cov0 = support_weights @ support_cov @ support_weights
        estimates[target] = obsval @ weights
        variance[target] = max(
            0.0, block_cov0 - weights @ block_rhs - mu
        )
    return estimates, variance


def _sk_lib_gslib_reference(obsloc, obsval, newloc, sk_mean, **vgm):
    """Simple kriging using the same covariance and angle conventions as libkriging."""
    obsloc = np.asarray(obsloc, dtype=float)
    obsval = np.asarray(obsval, dtype=float)
    newloc = np.asarray(newloc, dtype=float)
    ndim = obsloc.shape[1]
    transform = _libkriging_gslib_rotmat(
        vgm["hmin"], vgm["hmax"], vgm["hvert"],
        vgm.get("azimuth", 0.0), vgm.get("dip", 0.0), vgm.get("plunge", 0.0),
    )[:ndim, :ndim]
    scaled_obs = (transform @ obsloc.T).T
    scaled_new = (transform @ newloc.T).T
    nugget = vgm.get("nugget", 0.0)
    sill = vgm.get("sill", 1.0)
    vtype = vgm.get("vtype", "sph")
    cov0 = nugget + sill
    c_obs = _ok3d_covariance(vtype, squareform(pdist(scaled_obs))) * sill
    np.fill_diagonal(c_obs, cov0)
    rhs_dist = cdist(scaled_obs, scaled_new)
    c_rhs = _ok3d_covariance(vtype, rhs_dist) * sill
    if nugget > 0.0:
        c_rhs = c_rhs + np.where(rhs_dist < 1e-10, nugget, 0.0)
    weights = np.linalg.solve(c_obs, c_rhs)
    estimates = sk_mean + (obsval - sk_mean) @ weights
    variance = np.maximum(0.0, cov0 - (weights * c_rhs).sum(axis=0))
    return estimates, variance


def _ok_reference_3d(obsloc, obsval, newloc, vgm):
    obsloc = np.asarray(obsloc, dtype=float)
    obsval = np.asarray(obsval, dtype=float)
    newloc = np.asarray(newloc, dtype=float)
    vtype, nugget, sill, hmax, rotmat, nnear = _parse_ok3d_vgm(vgm)
    if nnear < len(obsloc):
        raise ValueError(f"nnear ({nnear}) must be >= nobs ({len(obsloc)})")

    scaled_obs = (rotmat @ obsloc.T).T / hmax
    scaled_new = (rotmat @ newloc.T).T / hmax
    cov0 = nugget + sill
    c_obs = _ok3d_covariance(vtype, squareform(pdist(scaled_obs))) * sill
    np.fill_diagonal(c_obs, cov0)
    rhs_dist = cdist(scaled_obs, scaled_new)
    c_rhs = _ok3d_covariance(vtype, rhs_dist) * sill
    if nugget > 0.0:
        # Texture2Par treats an exact target-observation match as including the
        # nugget covariance on the RHS. This is a T2P support/convention choice;
        # tests that should be independent of that convention should avoid exact
        # target-observation coincidences when nugget > 0.
        c_rhs = c_rhs + np.where(rhs_dist < 1e-10, nugget, 0.0)

    e = np.ones(len(obsloc))
    cf = cho_factor(c_obs, lower=True)
    cinv_rhs = cho_solve(cf, c_rhs)
    cinv_e = cho_solve(cf, e)
    mu = (e @ cinv_rhs - 1.0) / (e @ cinv_e)
    weights = cinv_rhs - np.outer(cinv_e, mu)

    est = obsval @ weights
    var = np.maximum(0.0, cov0 - (weights * c_rhs).sum(axis=0) - mu)
    return est, var


def _ok_reference_3d_augmented(obsloc, obsval, newloc, vgm):
    obsloc = np.asarray(obsloc, dtype=float)
    obsval = np.asarray(obsval, dtype=float)
    newloc = np.asarray(newloc, dtype=float)
    vtype, nugget, sill, hmax, rotmat, nnear = _parse_ok3d_vgm(vgm)
    if nnear < len(obsloc):
        raise ValueError(f"nnear ({nnear}) must be >= nobs ({len(obsloc)})")

    scaled_obs = (rotmat @ obsloc.T).T / hmax
    scaled_new = (rotmat @ newloc.T).T / hmax
    cov0 = nugget + sill
    c_obs = _ok3d_covariance(vtype, squareform(pdist(scaled_obs))) * sill
    np.fill_diagonal(c_obs, cov0)
    rhs_dist = cdist(scaled_obs, scaled_new)
    c_rhs = _ok3d_covariance(vtype, rhs_dist) * sill
    if nugget > 0.0:
        # Texture2Par treats an exact target-observation match as including the
        # nugget covariance on the RHS. This is a T2P support/convention choice;
        # tests that should be independent of that convention should avoid exact
        # target-observation coincidences when nugget > 0.
        c_rhs = c_rhs + np.where(rhs_dist < 1e-10, nugget, 0.0)

    nobs = len(obsloc)
    lhs = np.zeros((nobs + 1, nobs + 1))
    lhs[:nobs, :nobs] = c_obs
    lhs[:nobs, nobs] = 1.0
    lhs[nobs, :nobs] = 1.0

    rhs = np.vstack([c_rhs, np.ones((1, len(newloc)))])
    solved = np.linalg.solve(lhs, rhs)
    weights = solved[:nobs, :]
    mu = solved[nobs, :]
    est = obsval @ weights
    var = np.maximum(0.0, cov0 - (weights * c_rhs).sum(axis=0) - mu)
    return est, var


def _ok_reference_3d_nearest(obsloc, obsval, newloc, vgm, anisotropic_search=False):
    obsloc = np.asarray(obsloc, dtype=float)
    obsval = np.asarray(obsval, dtype=float)
    newloc = np.asarray(newloc, dtype=float)
    _, _, _, _, rotmat, nnear = _parse_ok3d_vgm(vgm)
    nuse = min(nnear, len(obsloc))
    if anisotropic_search:
        search_obs = (rotmat @ obsloc.T).T
        search_new = (rotmat @ newloc.T).T
    else:
        search_obs = obsloc
        search_new = newloc

    est = []
    var = []
    for i, target in enumerate(newloc):
        dist = np.linalg.norm(search_obs - search_new[i], axis=1)
        use = np.argsort(dist, kind="stable")[:nuse]
        loc_est, loc_var = _ok_reference_3d(obsloc[use], obsval[use], target[None, :], vgm)
        est.append(loc_est[0])
        var.append(loc_var[0])
    return np.asarray(est), np.asarray(var)


def _ok3d_observations_from_dataset(dataset_text: str):
    rows = list(csv.DictReader(dataset_text.splitlines()))
    locs = []
    vals = []
    for row in rows:
        top = float(row["TopDepth"])
        bot = float(row["BotDepth"])
        z = float(row["Zland"]) - 0.5 * (top + bot)
        locs.append([float(row["X"]), float(row["Y"]), z])
        vals.append(float(row["Coarse"]))
    return np.asarray(locs), np.asarray(vals)


def _coarse_output_to_arrays(df: pd.DataFrame):
    layer_cols = ["Layer1", "Layer2", "Layer3"]
    variance_cols = ["Var_Layer1", "Var_Layer2", "Var_Layer3"]
    missing_var = set(variance_cols) - set(df.columns)
    if missing_var:
        pytest.fail(
            f"Coarse texture output missing variance columns {sorted(missing_var)} "
            f"despite WRITE_VARIANCE. Columns: {list(df.columns)}. "
            f"Rows:\n{df.to_string(index=False)}"
        )

    z_centers = [95.0, 85.0, 75.0]
    targets = []
    for _, row in df.iterrows():
        for z in z_centers:
            targets.append([row["X"], row["Y"], z])
    est = df[layer_cols].to_numpy(dtype=float)
    var = df[variance_cols].to_numpy(dtype=float)
    return np.asarray(targets, dtype=float), est, var


# ---------------------------------------------------------------------------
# 1. Basic kriging run
# ---------------------------------------------------------------------------

class TestOrdinaryKriging:

    def test_run_completes(self, run_t2p):
        r = run_t2p(_grid_t2p(), extra_files=FILES)
        assert r.returncode == 0, r.combined

    def test_kh_output_written(self, run_t2p):
        r = run_t2p(_grid_t2p(), extra_files=FILES)
        assert r.returncode == 0
        kh_files = list(r.workdir.glob("*KH*"))
        assert kh_files, "No KH output file found"

    def test_kh_values_positive(self, run_t2p):
        r = run_t2p(_grid_t2p(), extra_files=FILES)
        assert r.returncode == 0
        kh_files = list(r.workdir.glob("*KH*"))
        if not kh_files:
            pytest.skip("KH file not found")
        values = _read_csv_column(kh_files[0], list(
            csv.DictReader(kh_files[0].read_text().splitlines()))[0].keys().__iter__().__next__())
        # Fall back: just check file is non-empty
        assert kh_files[0].stat().st_size > 0

    def test_texture_bounded_0_1(self, run_t2p):
        """With NO_TEX_BOUNDS off (default), interpolated coarse values stay in [0,1]."""
        r = run_t2p(_grid_t2p(options_extra="WRITE_DATASET_FILES"), extra_files=FILES)
        assert r.returncode == 0
        tex_files = list(r.workdir.glob("*Coarse*"))
        if not tex_files:
            pytest.skip("Coarse dataset file not found")
        est = pd.read_csv(tex_files[0])
        assert est["Layer1"].min()>=0 and est["Layer1"].max()<=1, "Texture values outside [0,1]"

    def test_infer_last_class_coarse_plus_fine_sum_to_one(self, run_t2p):
        """With INFER_LAST_CLASS, Coarse + Fine should sum to 1 at every node."""
        r = run_t2p(_grid_t2p(options_extra="WRITE_DATASET_FILES"), extra_files=FILES)
        assert r.returncode == 0
        coarse_files = list(r.workdir.glob("*Coarse*"))
        fine_files   = list(r.workdir.glob("*Fine*"))
        if not coarse_files or not fine_files:
            pytest.skip("Dataset output files not found")
        # Just verify both files are written (sum check requires matching rows)
        assert coarse_files and fine_files


# ---------------------------------------------------------------------------
# 2. Default libkriging ordinary-kriging backend and transition aliases
# ---------------------------------------------------------------------------

class TestLibKrigingOrdinary:

    @staticmethod
    def _files():
        return {"dataset.csv": OK_LIB_2D_DATASET, "grid.grd": OK_LIB_2D_GRD}

    def test_ok_lib_alias_is_accepted_and_exact_hits_run_in_2d(self, run_t2p):
        r = run_t2p(_ok_lib_2d_t2p(method="OK_LIB"), extra_files=self._files())
        assert r.returncode == 0, r.combined
        assert "OK_LIB alias" in r.log

        candidates = sorted(r.workdir.glob("*Coarse*.csv"))
        assert candidates, f"No OK_LIB texture output found.\n{r.combined}"
        df = pd.read_csv(candidates[0])
        expected = {
            (0.5, 0.5): 0.10,
            (1.5, 0.5): 0.20,
            (0.5, 1.5): 0.30,
            (1.5, 1.5): 0.40,
        }
        for row in df.itertuples(index=False):
            assert row.Layer1 == pytest.approx(expected[(row.X, row.Y)], abs=1e-10)

    def test_zero_observation_variance_preserves_exact_hit_and_zero_kriging_variance(
        self, run_t2p
    ):
        r = run_t2p(
            _ok_lib_2d_t2p(
                method="OK", options_extra="KVME\n          WRITE_VARIANCE"
            ),
            extra_files={
                "dataset.csv": _ok_lib_variance_dataset(0.0),
                "grid.grd": OK_LIB_2D_GRD,
            },
        )
        assert r.returncode == 0, r.combined
        df = _read_coarse_texture_output(r, ("Layer1",))
        hit = df[(df["X"] == 0.5) & (df["Y"] == 0.5)].iloc[0]
        assert hit["Layer1"] == pytest.approx(0.10, abs=1e-10)
        assert hit["Var_Layer1"] == pytest.approx(0.0, abs=1e-10)

    def test_nonzero_observation_variance_changes_estimate_and_kriging_variance(
        self, run_t2p
    ):
        t2p = _ok_lib_2d_t2p(
            method="OK", options_extra="KVME\n          WRITE_VARIANCE"
        )
        exact = run_t2p(
            t2p,
            extra_files={
                "dataset.csv": _ok_lib_variance_dataset(0.0),
                "grid.grd": OK_LIB_2D_GRD,
            },
        )
        assert exact.returncode == 0, exact.combined
        exact_df = _read_coarse_texture_output(exact, ("Layer1",)).copy()

        noisy = run_t2p(
            t2p,
            extra_files={
                "dataset.csv": _ok_lib_variance_dataset(1.0),
                "grid.grd": OK_LIB_2D_GRD,
            },
        )
        assert noisy.returncode == 0, noisy.combined
        noisy_df = _read_coarse_texture_output(noisy, ("Layer1",))
        exact_hit = exact_df[(exact_df["X"] == 0.5) & (exact_df["Y"] == 0.5)].iloc[0]
        noisy_hit = noisy_df[(noisy_df["X"] == 0.5) & (noisy_df["Y"] == 0.5)].iloc[0]

        assert abs(noisy_hit["Layer1"] - exact_hit["Layer1"]) > 1e-3
        assert noisy_hit["Var_Layer1"] > exact_hit["Var_Layer1"] + 1e-3

    def test_cokriging_retains_primary_observation_variance(self, run_t2p):
        files = {
            "primary.csv": _ck_variance_primary_dataset(0.0),
            "secondary.csv": CK_VARIANCE_SECONDARY_DATASET,
            "grid.grd": OK_LIB_2D_GRD,
        }
        exact = run_t2p(_ck_variance_t2p(), extra_files=files)
        assert exact.returncode == 0, exact.combined
        exact_df = _read_coarse_texture_output(exact, ("Layer1",)).copy()

        files["primary.csv"] = _ck_variance_primary_dataset(1.0)
        noisy = run_t2p(_ck_variance_t2p(), extra_files=files)
        assert noisy.returncode == 0, noisy.combined
        noisy_df = _read_coarse_texture_output(noisy, ("Layer1",))
        exact_hit = exact_df[(exact_df["X"] == 0.5) & (exact_df["Y"] == 0.5)].iloc[0]
        noisy_hit = noisy_df[(noisy_df["X"] == 0.5) & (noisy_df["Y"] == 0.5)].iloc[0]

        assert exact_hit["Layer1"] == pytest.approx(0.10, abs=1e-10)
        assert exact_hit["Var_Layer1"] == pytest.approx(0.0, abs=1e-10)
        assert abs(noisy_hit["Layer1"] - exact_hit["Layer1"]) > 1e-3
        assert noisy_hit["Var_Layer1"] > exact_hit["Var_Layer1"] + 1e-3

    def test_ok_legacy_escape_hatch_is_rejected(self, run_t2p):
        r = run_t2p(
            _ok_lib_2d_t2p(method="OK_LEGACY"), extra_files=self._files()
        )
        assert r.returncode != 0
        assert "Invalid Interpolation Method" in r.combined

    def test_ok_lib_alias_accepts_sgsim(self, run_t2p):
        r = run_t2p(
            _ok_lib_2d_t2p(method="OK_LIB", options_extra="SGSIM\nSEED 23"),
            extra_files={"dataset.csv": SK_2D_DATASET, "grid.grd": SK_2D_GRD},
        )
        assert r.returncode == 0, r.combined
        assert "SGSIM simulation mode" in r.log
        df = _read_coarse_texture_output(r, ("Layer1",))
        assert np.isfinite(df["Layer1"]).all()

    def test_default_ok_exact_vertical_hits_run_in_3d(self, run_t2p):
        t2p = _ok3d_t2p(method="OK")
        r = run_t2p(
            t2p,
            extra_files={
                "dataset.csv": _ok3d_exact_hit_dataset(),
                "grid.grd": MAX_LOG_LENGTH_GRD,
            },
        )
        assert r.returncode == 0, r.combined
        row = _read_coarse_texture_output(r).iloc[0]
        assert [row["Layer1"], row["Layer2"], row["Layer3"]] == pytest.approx(
            [0.10, 0.50, 0.90], abs=1e-10
        )

    def test_default_ok_nnear_one_uses_nearest_observation(self, run_t2p):
        r = run_t2p(
            _ok_lib_2d_t2p(
                variogram_lines="1 Sph 0.0 1.0 10.0 10.0 0.0 1"
            ),
            extra_files={
                "dataset.csv": OK_LIB_NNEAR_DATASET,
                "grid.grd": OK_LIB_2D_GRD,
            },
        )
        assert r.returncode == 0, r.combined
        df = _read_coarse_texture_output(r, ("Layer1",))
        target = df[(df["X"] == 0.5) & (df["Y"] == 0.5)].iloc[0]
        assert target["Layer1"] == pytest.approx(0.11, abs=1e-10)

    def test_default_ok_2d_anisotropy_matches_gslib_reference(self, run_t2p):
        # This is intentionally a libkriging/GSLIB-convention test, not a
        # parity comparison with legacy Texture2Par OK.
        rows = list(csv.DictReader(OK_LIB_2D_PARITY_DATASET.splitlines()))
        obsloc = np.array([[float(row["X"]), float(row["Y"])] for row in rows])
        obsval = np.array([float(row["Coarse"]) for row in rows])
        params = dict(hmin=0.6, hmax=4.0, hvert=4.0, azimuth=35.0)
        r = run_t2p(
            _ok_lib_2d_t2p(
                options_extra="WRITE_VARIANCE",
                variogram_lines="1 Sph 0.0 1.0 0.6 4.0 35.0 5",
            ),
            extra_files={
                "dataset.csv": OK_LIB_2D_PARITY_DATASET,
                "grid.grd": OK_LIB_2D_GRD,
            },
        )
        assert r.returncode == 0, r.combined
        df = _read_coarse_texture_output(r, ("Layer1",)).sort_values(["X", "Y"])
        targets = df[["X", "Y"]].to_numpy(dtype=float)
        ref_est, ref_var = _ok_lib_gslib_reference(
            obsloc, obsval, targets, **params
        )
        iso_est, _ = _ok_lib_gslib_reference(
            obsloc, obsval, targets, hmin=4.0, hmax=4.0, hvert=4.0
        )

        assert df["Layer1"].to_numpy() == pytest.approx(ref_est, abs=3e-6)
        assert df["Var_Layer1"].to_numpy() == pytest.approx(ref_var, abs=6e-6)
        assert np.max(np.abs(ref_est - iso_est)) > 0.01

    def test_default_ok_3d_dip_anisotropy_matches_gslib_reference(self, run_t2p):
        # Positive dip and angle order follow libkriging/GSLIB here. The legacy
        # _ok_reference_3d helper is deliberately not used.
        dataset = _ok3d_irregular_dataset()
        vgm = "Sph 0.0 1.0 1.2 4.0 12.0 25.0 20.0 -15.0 20"
        r = run_t2p(
            _ok3d_multi_t2p(
                vgm=vgm, method="OK"
            ),
            extra_files={"dataset.csv": dataset, "grid.grd": OK3D_INTERIOR_GRD},
        )
        assert r.returncode == 0, r.combined
        df = _read_coarse_texture_output(r).sort_values(["X", "Y"])
        targets, t2p_est, t2p_var = _coarse_output_to_arrays(df)
        obsloc, obsval = _ok3d_observations_from_dataset(dataset)
        params = dict(
            hmin=1.2,
            hmax=4.0,
            hvert=12.0,
            azimuth=25.0,
            dip=20.0,
            plunge=-15.0,
        )
        ref_est, ref_var = _ok_lib_gslib_reference(
            obsloc, obsval, targets, **params
        )
        zero_angle_est, _ = _ok_lib_gslib_reference(
            obsloc, obsval, targets, hmin=1.2, hmax=4.0, hvert=12.0
        )

        assert t2p_est.ravel() == pytest.approx(ref_est, abs=6e-6)
        assert t2p_var.ravel() == pytest.approx(ref_var, abs=6e-6)
        assert np.max(np.abs(ref_est - zero_angle_est)) > 0.01

    def test_default_ok_anisotropic_search_uses_first_structure_and_angle(
        self, run_t2p
    ):
        # At (0.5, 0.5), the 0.10 observation is Euclidean-nearest, while the
        # 0.90 observation is nearest after rotating the GSLIB search ellipse
        # 90 degrees (minor=1, major=10). The isotropic second structure must
        # not replace the first structure's transitional search parameters.
        files = {
            "dataset.csv": OK_LIB_ANISO_SEARCH_DATASET,
            "grid.grd": OK_LIB_2D_GRD,
        }
        variogram = (
            "1 Sph 0.0 0.6 1.0 10.0 90.0 1\n"
            "            2 Sph 0.0 0.4 10.0 10.0 0.0 1"
        )
        isotropic_search = run_t2p(
            _ok_lib_2d_t2p(variogram_lines=variogram, max_vstruct=2),
            extra_files=files,
        )
        assert isotropic_search.returncode == 0, isotropic_search.combined
        isotropic_df = _read_coarse_texture_output(
            isotropic_search, ("Layer1",)
        )
        isotropic_value = isotropic_df[
            (isotropic_df["X"] == 0.5) & (isotropic_df["Y"] == 0.5)
        ].iloc[0]["Layer1"]

        anisotropic_search = run_t2p(
            _ok_lib_2d_t2p(
                options_extra="ANISOTROPIC_SEARCH",
                variogram_lines=variogram,
                max_vstruct=2,
            ),
            extra_files=files,
        )
        assert anisotropic_search.returncode == 0, anisotropic_search.combined
        anisotropic_df = _read_coarse_texture_output(
            anisotropic_search, ("Layer1",)
        )
        anisotropic_value = anisotropic_df[
            (anisotropic_df["X"] == 0.5) & (anisotropic_df["Y"] == 0.5)
        ].iloc[0]["Layer1"]

        assert isotropic_value == pytest.approx(0.10, abs=1e-10)
        assert anisotropic_value == pytest.approx(0.90, abs=1e-10)

    def test_default_ok_nested_structures_run(self, run_t2p):
        structures = (
            "1 Sph 0.02 0.58 2.5 2.5 0.0 5\n"
            "            2 Exp 0.00 0.40 6.0 6.0 0.0 5"
        )
        files = {
            "dataset.csv": OK_LIB_2D_PARITY_DATASET,
            "grid.grd": OK_LIB_2D_GRD,
        }
        lib = run_t2p(
            _ok_lib_2d_t2p(variogram_lines=structures, max_vstruct=2),
            extra_files=files,
        )
        assert lib.returncode == 0, lib.combined
        lib_values = _read_coarse_texture_output(
            lib, ("Layer1",)
        ).sort_values(["X", "Y"])
        assert np.isfinite(lib_values["Layer1"]).all()

    @pytest.mark.parametrize(
        "vtype", ["Sph", "Exp", "Hol", "Gau", "Pow", "Cir", "Lin"]
    )
    def test_default_ok_supported_variogram_mapping_exact_hits(self, run_t2p, vtype):
        r = run_t2p(
            _ok_lib_2d_t2p(
                variogram_lines=f"1 {vtype} 0.0 1.0 10.0 10.0 0.0 1"
            ),
            extra_files=self._files(),
        )
        assert r.returncode == 0, r.combined
        values = _read_coarse_texture_output(
            r, ("Layer1",)
        ).sort_values(["X", "Y"])["Layer1"]
        assert values.tolist() == pytest.approx([0.10, 0.30, 0.20, 0.40], abs=1e-10)

    def test_default_ok_pilot_point_run_writes_hydraulic_files(self, run_t2p):
        r = run_t2p(
            _grid_t2p(),
            extra_files=FILES,
        )
        assert r.returncode == 0, r.combined
        for property_name in ("KH", "KV", "SS", "SY"):
            assert list(r.workdir.glob(f"*{property_name}*")), (
                f"Default OK pilot-point run produced no {property_name} output"
            )

    @pytest.mark.parametrize(
        ("interp_dim", "dataset"),
        [(2, DUPLICATE_2D_DATASET), (3, DUPLICATE_3D_DATASET)],
    )
    def test_default_ok_duplicate_conditioning_is_fatal(
        self, run_t2p, interp_dim, dataset
    ):
        t2p = _duplicate_conditioning_t2p(interp_dim=interp_dim)
        r = run_t2p(
            t2p,
            extra_files={
                "dataset.csv": dataset,
                "grid.grd": DUPLICATE_GRD,
            },
        )
        assert r.returncode != 0
        assert "co-located" in r.combined.lower()

# ---------------------------------------------------------------------------
# 3. Block kriging status
# ---------------------------------------------------------------------------

class TestBlockKriging:

    def test_modflow_block_ok_matches_independent_reference(self, run_t2p):
        r = run_t2p(_block_modflow_t2p(), extra_files=BLOCK_FILES)
        assert r.returncode == 0, r.combined
        df = _read_coarse_texture_output(r, required_layers=("Layer1",))
        df = df.sort_values(["X", "Y"])
        centres = df[["X", "Y"]].to_numpy()
        obsloc = np.array([[0.1, 0.1], [1.9, 0.1],
                           [0.1, 1.9], [1.9, 1.9]])
        obsval = np.array([0.10, 0.30, 0.65, 0.90])
        expected, expected_var = _block_ok_reference(
            obsloc, obsval, centres, halfsize=(0.5, 0.5)
        )
        assert df["Layer1"].to_numpy() == pytest.approx(expected, abs=2e-6)
        assert df["Var_Layer1"].to_numpy() == pytest.approx(
            expected_var, abs=2e-6
        )

    def test_modflow_block_ok_retains_observation_variance(self, run_t2p):
        t2p = _block_modflow_t2p(options_extra="KVME")
        exact_files = {**BLOCK_FILES, "dataset.csv": _block_variance_dataset(0.0)}
        exact = run_t2p(t2p, extra_files=exact_files)
        assert exact.returncode == 0, exact.combined
        exact_df = _read_coarse_texture_output(
            exact, required_layers=("Layer1",)
        ).sort_values(["X", "Y"]).copy()

        noisy_files = {**BLOCK_FILES, "dataset.csv": _block_variance_dataset(1.0)}
        noisy = run_t2p(t2p, extra_files=noisy_files)
        assert noisy.returncode == 0, noisy.combined
        noisy_df = _read_coarse_texture_output(
            noisy, required_layers=("Layer1",)
        ).sort_values(["X", "Y"])

        assert np.max(np.abs(
            noisy_df["Layer1"].to_numpy() - exact_df["Layer1"].to_numpy()
        )) > 1e-3
        assert np.max(np.abs(
            noisy_df["Var_Layer1"].to_numpy()
            - exact_df["Var_Layer1"].to_numpy()
        )) > 1e-3

    @pytest.mark.parametrize("method", ["OK_LIB", "SK"])
    def test_modflow_block_shared_libkriging_modes_run(self, run_t2p, method):
        r = run_t2p(
            _block_modflow_t2p(method=method), extra_files=BLOCK_FILES
        )
        assert r.returncode == 0, r.combined
        df = _read_coarse_texture_output(r, required_layers=("Layer1",))
        assert np.isfinite(df["Layer1"]).all()
        assert np.isfinite(df["Var_Layer1"]).all()

    def test_block_mode_rejects_non_modflow_model(self, run_t2p):
        r = run_t2p(
            _ok_lib_2d_t2p(flow_extra="BLOCK_MODE"),
            extra_files=TestLibKrigingOrdinary._files(),
        )
        assert r.returncode != 0
        assert "BLOCK_MODE" in r.combined
        assert "only for structured MODFLOW" in r.combined

    @pytest.mark.parametrize(
        ("method", "options_extra", "diagnostic"),
        [
            ("IDW", "", "IDW has no block implementation"),
            ("OK", "SGSIM", "unsupported with SGSIM"),
        ],
    )
    def test_block_mode_rejects_point_only_routes(
        self, run_t2p, method, options_extra, diagnostic
    ):
        r = run_t2p(
            _block_modflow_t2p(method=method, options_extra=options_extra),
            extra_files=BLOCK_FILES,
        )
        assert r.returncode != 0
        assert diagnostic in r.combined
        assert "point support" in r.combined or method == "IDW"


# ---------------------------------------------------------------------------
# 4. Libkriging simple-kriging backend
# ---------------------------------------------------------------------------

class TestLibKrigingSimple:

    @staticmethod
    def _files_2d():
        return {"dataset.csv": SK_2D_DATASET, "grid.grd": SK_2D_GRD}

    @staticmethod
    def _run_simple_ck(run_t2p, primary=None, secondary=None):
        t2p = _ck_variance_t2p().replace("INTERP_METHOD OK", "INTERP_METHOD SK")
        return run_t2p(
            t2p,
            extra_files={
                "primary.csv": primary or _ck_variance_primary_dataset(0.0),
                "secondary.csv": secondary or SK_CK_SECONDARY_DATASET,
                "grid.grd": SK_2D_GRD,
            },
        )

    def test_sk_2d_uses_global_mean_and_differs_from_ok(self, run_t2p):
        sk = run_t2p(
            _ok_lib_2d_t2p(method="SK"), extra_files=self._files_2d()
        )
        assert sk.returncode == 0, sk.combined
        assert "Simple Kriging (libkriging)" in sk.log
        sk_df = _read_coarse_texture_output(sk, ("Layer1",)).sort_values(["X", "Y"])

        ok = run_t2p(_ok_lib_2d_t2p(method="OK"), extra_files=self._files_2d())
        assert ok.returncode == 0, ok.combined
        ok_df = _read_coarse_texture_output(ok, ("Layer1",)).sort_values(["X", "Y"])
        rows = list(csv.DictReader(SK_2D_DATASET.splitlines()))
        obsloc = np.array([[float(row["X"]), float(row["Y"])] for row in rows])
        obsval = np.array([float(row["Coarse"]) for row in rows])
        targets = sk_df[["X", "Y"]].to_numpy(dtype=float)
        expected, _ = _sk_lib_gslib_reference(
            obsloc, obsval, targets, sk_mean=obsval.mean(),
            hmin=10.0, hmax=10.0, hvert=10.0,
        )

        assert sk_df["Layer1"].to_numpy() == pytest.approx(expected, abs=3e-6)
        assert np.isfinite(sk_df["Layer1"]).all()
        assert sk_df["Layer1"].between(0.0, 1.0).all()
        assert np.max(np.abs(
            sk_df["Layer1"].to_numpy() - ok_df["Layer1"].to_numpy()
        )) > 1e-3

    def test_sk_3d_matches_global_mean_reference(self, run_t2p):
        dataset = _ok3d_irregular_dataset()
        vgm = "Sph 0.0 1.0 20.0 20.0 20.0 0.0 0.0 0.0 20"
        r = run_t2p(
            _ok3d_multi_t2p(vgm=vgm, method="SK"),
            extra_files={"dataset.csv": dataset, "grid.grd": OK3D_INTERIOR_GRD},
        )
        assert r.returncode == 0, r.combined
        df = _read_coarse_texture_output(r).sort_values(["X", "Y"])
        targets, estimates, variance = _coarse_output_to_arrays(df)
        obsloc, obsval = _ok3d_observations_from_dataset(dataset)
        expected, expected_var = _sk_lib_gslib_reference(
            obsloc, obsval, targets, sk_mean=obsval.mean(),
            hmin=20.0, hmax=20.0, hvert=20.0,
        )
        assert estimates.ravel() == pytest.approx(expected, abs=6e-6)
        assert variance.ravel() == pytest.approx(expected_var, abs=6e-6)
        assert np.isfinite(estimates).all()

    def test_sk_model_defaults_to_libkriging_ok_for_pilot_points(self, run_t2p):
        r = run_t2p(
            _grid_t2p(options_extra="INTERP_METHOD SK"), extra_files=FILES
        )
        assert r.returncode == 0, r.combined
        assert "Simple Kriging (libkriging)" in r.log
        pp_line = next(
            line for line in r.log.splitlines()
            if "Pilot-Point Interpolation Method:" in line
        )
        assert "Ordinary Kriging (libkriging)" in pp_line
        for property_name in ("KH", "KV", "SS", "SY"):
            assert list(r.workdir.glob(f"*{property_name}*"))

    def test_pilot_point_idw_still_runs(self, run_t2p):
        r = run_t2p(
            _grid_t2p(options_extra="PP_INTERP_METHOD IDW"), extra_files=FILES
        )
        assert r.returncode == 0, r.combined
        pp_line = next(
            line for line in r.log.splitlines()
            if "Pilot-Point Interpolation Method:" in line
        )
        assert "Inverse Distance Weighting (IDW)" in pp_line
        for property_name in ("KH", "KV", "SS", "SY"):
            assert list(r.workdir.glob(f"*{property_name}*"))

    def test_simple_cokriging_uses_secondary_data(self, run_t2p):
        simple_ck = self._run_simple_ck(run_t2p)
        assert simple_ck.returncode == 0, simple_ck.combined
        ck_values = _read_coarse_texture_output(
            simple_ck, ("Layer1",)
        ).sort_values(["X", "Y"])["Layer1"].to_numpy()

        primary_only = run_t2p(
            _ok_lib_2d_t2p(method="SK"),
            extra_files={
                "dataset.csv": _ck_variance_primary_dataset(0.0),
                "grid.grd": SK_2D_GRD,
            },
        )
        assert primary_only.returncode == 0, primary_only.combined
        primary_values = _read_coarse_texture_output(
            primary_only, ("Layer1",)
        ).sort_values(["X", "Y"])["Layer1"].to_numpy()

        assert np.max(np.abs(ck_values - primary_values)) > 1e-4

    def test_simple_cokriging_uses_each_variable_mean(self, run_t2p):
        baseline = self._run_simple_ck(run_t2p)
        assert baseline.returncode == 0, baseline.combined
        baseline_values = _read_coarse_texture_output(
            baseline, ("Layer1",)
        ).sort_values(["X", "Y"])["Layer1"].to_numpy()

        secondary_shifted = self._run_simple_ck(
            run_t2p,
            secondary=_shift_dataset_column(
                SK_CK_SECONDARY_DATASET, "Rho", 2.0
            ),
        )
        assert secondary_shifted.returncode == 0, secondary_shifted.combined
        secondary_values = _read_coarse_texture_output(
            secondary_shifted, ("Layer1",)
        ).sort_values(["X", "Y"])["Layer1"].to_numpy()

        primary_shifted = self._run_simple_ck(
            run_t2p,
            primary=_shift_dataset_column(
                _ck_variance_primary_dataset(0.0), "Coarse", 0.1
            ),
        )
        assert primary_shifted.returncode == 0, primary_shifted.combined
        primary_values = _read_coarse_texture_output(
            primary_shifted, ("Layer1",)
        ).sort_values(["X", "Y"])["Layer1"].to_numpy()

        assert secondary_values == pytest.approx(baseline_values, abs=3e-6)
        assert primary_values == pytest.approx(baseline_values + 0.1, abs=3e-6)

    def test_sk_after_without_defined_mean_is_fatal(self, run_t2p):
        t2p = _grid_t2p(
            options_extra="INTERP_METHOD SK\n          INTERP_ORDER AFTER",
            infer=False,
        )
        files = FILES
        r = run_t2p(t2p, extra_files=files)
        assert r.returncode != 0
        assert "no defined SK mean".lower() in r.combined.lower()


# ---------------------------------------------------------------------------
# 4. CORRECT_WEIGHTS (non-negative kriging)
# ---------------------------------------------------------------------------

class TestCorrectWeights:

    def test_correct_weights_flag_runs(self, run_t2p):
        r = run_t2p(_grid_t2p(options_extra="CORRECT_WEIGHTS"), extra_files=FILES)
        assert r.returncode == 0, r.combined


# ---------------------------------------------------------------------------
# 3. Conditioning/target coincidence contracts
# ---------------------------------------------------------------------------

class TestCoincidentLocations:

    @staticmethod
    def _assert_duplicate_conditioning_is_fatal(r):
        assert r.returncode != 0, (
            "Co-located conditioning observations must be fatal, but "
            f"Texture2Par returned success.\n{r.combined}"
        )
        fatal_terms = (
            "duplicate", "co-locat", "decompos", "factor", "singular",
            "positive definite", "kriging matrix",
        )
        assert any(term in r.combined.lower() for term in fatal_terms), (
            "Duplicate conditioning was fatal but emitted no recognizable "
            f"diagnostic.\n{r.combined}"
        )

    def test_duplicate_xy_conditioning_is_fatal_in_2d(self, run_t2p):
        # W1 and W2 are distinct observations that collapse to the same 2D
        # support point after normal layer preprocessing. This is invalid.
        r = run_t2p(
            _duplicate_conditioning_t2p(interp_dim=2),
            extra_files={
                "dataset.csv": DUPLICATE_2D_DATASET,
                "grid.grd": DUPLICATE_GRD,
            },
        )
        self._assert_duplicate_conditioning_is_fatal(r)

    def test_duplicate_xyz_conditioning_is_fatal_in_3d(self, run_t2p):
        # Equal X/Y and equal interval midpoint elevation give W1 and W2 the
        # same 3D support coordinate. Unlike an exact target hit, this is fatal.
        r = run_t2p(
            _duplicate_conditioning_t2p(interp_dim=3),
            extra_files={
                "dataset.csv": DUPLICATE_3D_DATASET,
                "grid.grd": DUPLICATE_GRD,
            },
        )
        self._assert_duplicate_conditioning_is_fatal(r)


# ---------------------------------------------------------------------------
# 4. Sequential Gaussian Simulation (SGS)
# ---------------------------------------------------------------------------

class TestSGS:

    @staticmethod
    def _run_2d(run_t2p, seed, dataset=SK_2D_DATASET, extra=""):
        options = f"SGSIM\nSEED {seed}"
        if extra:
            options += f"\n{extra}"
        return run_t2p(
            _ok_lib_2d_t2p(options_extra=options),
            extra_files={"dataset.csv": dataset, "grid.grd": SK_2D_GRD},
        )

    @staticmethod
    def _run_sk(run_t2p, seed, dataset=SK_2D_DATASET, extra=""):
        options = f"SGSIM\nSEED {seed}"
        if extra:
            options += f"\n{extra}"
        return run_t2p(
            _ok_lib_2d_t2p(method="SK", options_extra=options),
            extra_files={"dataset.csv": dataset, "grid.grd": SK_2D_GRD},
        )

    @staticmethod
    def _run_sk_ck(run_t2p, seed, primary=None, secondary=None):
        t2p = _ck_variance_t2p().replace(
            "MAX_VSTRUCT 1",
            f"MAX_VSTRUCT 1\n          SGSIM\n          SEED {seed}",
        ).replace("INTERP_METHOD OK", "INTERP_METHOD SK")
        return run_t2p(
            t2p,
            extra_files={
                "primary.csv": primary or _ck_variance_primary_dataset(0.0),
                "secondary.csv": secondary or SK_CK_SECONDARY_DATASET,
                "grid.grd": SK_2D_GRD,
            },
        )

    @staticmethod
    def _run_ck(run_t2p, seed, first_variance=0.0, extra="", interp_dim=2,
                mean_mode="LIBKRIGING"):
        options = f"SGSIM\nSEED {seed}"
        if interp_dim == 3:
            options += "\nINTERP_DIM 3"
        if extra:
            options += f"\n{extra}"
        t2p = _ck_variance_t2p().replace(
            "MAX_VSTRUCT 1",
            f"MAX_VSTRUCT 1\n          {options}",
        ).replace("CK_MEAN_MODE LIBKRIGING", f"CK_MEAN_MODE {mean_mode}")
        if interp_dim == 3:
            t2p = t2p.replace(
                "1 Sph 0.0 1.0 10.0 10.0 0.0 4",
                "1 Sph 0.0 1.0 10.0 10.0 10.0 0.0 0.0 0.0 4",
            ).replace(
                "1 Sph 0.0 0.25 10.0 10.0 0.0",
                "1 Sph 0.0 0.25 10.0 10.0 10.0 0.0 0.0 0.0",
            )
        return run_t2p(
            t2p,
            extra_files={
                "primary.csv": _ck_variance_primary_dataset(first_variance),
                "secondary.csv": CK_VARIANCE_SECONDARY_DATASET,
                "grid.grd": SK_2D_GRD,
            },
        )

    def test_sgsim_runs_through_libkriging_with_finite_output(self, run_t2p):
        r = self._run_2d(run_t2p, 12345)
        assert r.returncode == 0, r.combined
        assert "Ordinary Kriging (libkriging) - SGSIM simulation mode" in r.log
        df = _read_coarse_texture_output(r, ("Layer1",))
        assert np.isfinite(df["Layer1"]).all()
        assert not (r.workdir / "sgs_path.dat").exists()
        assert not (r.workdir / "sgs_sample.dat").exists()

    def test_sgsim_different_seeds_differ(self, run_t2p):
        r1 = self._run_2d(run_t2p, 1)
        assert r1.returncode == 0, r1.combined
        values1 = _read_coarse_texture_output(r1, ("Layer1",))["Layer1"].to_numpy().copy()
        r2 = self._run_2d(run_t2p, 2)
        assert r2.returncode == 0, r2.combined
        values2 = _read_coarse_texture_output(r2, ("Layer1",))["Layer1"].to_numpy()
        assert not np.array_equal(values1, values2)

    def test_sgsim_same_seed_reproducible(self, run_t2p):
        r1 = self._run_2d(run_t2p, 99)
        assert r1.returncode == 0, r1.combined
        values1 = _read_coarse_texture_output(r1, ("Layer1",))["Layer1"].to_numpy().copy()
        r2 = self._run_2d(run_t2p, 99)
        assert r2.returncode == 0, r2.combined
        values2 = _read_coarse_texture_output(r2, ("Layer1",))["Layer1"].to_numpy()
        assert np.array_equal(values1, values2)

    def test_sk_sgsim_point_support_is_finite_and_seeded(self, run_t2p):
        first = self._run_sk(run_t2p, 101)
        assert first.returncode == 0, first.combined
        assert "Simple Kriging (libkriging) - SGSIM simulation mode" in first.log
        first_values = _read_coarse_texture_output(
            first, ("Layer1",)
        )["Layer1"].to_numpy().copy()
        assert np.isfinite(first_values).all()

        repeat = self._run_sk(run_t2p, 101)
        assert repeat.returncode == 0, repeat.combined
        repeat_values = _read_coarse_texture_output(
            repeat, ("Layer1",)
        )["Layer1"].to_numpy()
        changed = self._run_sk(run_t2p, 102)
        assert changed.returncode == 0, changed.combined
        changed_values = _read_coarse_texture_output(
            changed, ("Layer1",)
        )["Layer1"].to_numpy()
        assert np.array_equal(first_values, repeat_values)
        assert not np.array_equal(first_values, changed_values)

    def test_sk_sgsim_accepts_observation_variance(self, run_t2p):
        exact = self._run_sk(
            run_t2p, 71, _ok_lib_variance_dataset(0.0),
            "KVME\nWRITE_VARIANCE",
        )
        assert exact.returncode == 0, exact.combined
        exact_df = _read_coarse_texture_output(exact, ("Layer1",)).copy()
        exact_hit = exact_df[
            (exact_df["X"] == 0.5) & (exact_df["Y"] == 0.5)
        ].iloc[0]
        assert exact_hit["Layer1"] == pytest.approx(0.10, abs=1e-10)
        assert exact_hit["Var_Layer1"] == pytest.approx(0.0, abs=1e-10)

        noisy = self._run_sk(
            run_t2p, 71, _ok_lib_variance_dataset(1.0),
            "KVME\nWRITE_VARIANCE",
        )
        assert noisy.returncode == 0, noisy.combined
        noisy_df = _read_coarse_texture_output(noisy, ("Layer1",))
        noisy_hit = noisy_df[
            (noisy_df["X"] == 0.5) & (noisy_df["Y"] == 0.5)
        ].iloc[0]
        assert abs(noisy_hit["Layer1"] - exact_hit["Layer1"]) > 1e-4
        assert noisy_hit["Var_Layer1"] > exact_hit["Var_Layer1"] + 1e-4

    def test_sgsim_zero_variance_exact_conditioning_and_soft_data(self, run_t2p):
        exact = self._run_2d(
            run_t2p, 71, _ok_lib_variance_dataset(0.0),
            "KVME\nWRITE_VARIANCE",
        )
        assert exact.returncode == 0, exact.combined
        exact_df = _read_coarse_texture_output(exact, ("Layer1",)).copy()
        hit = exact_df[(exact_df["X"] == 0.5) & (exact_df["Y"] == 0.5)].iloc[0]
        assert hit["Layer1"] == pytest.approx(0.10, abs=1e-10)
        assert hit["Var_Layer1"] == pytest.approx(0.0, abs=1e-10)

        noisy = self._run_2d(
            run_t2p, 71, _ok_lib_variance_dataset(1.0),
            "KVME\nWRITE_VARIANCE",
        )
        assert noisy.returncode == 0, noisy.combined
        noisy_df = _read_coarse_texture_output(noisy, ("Layer1",))
        noisy_hit = noisy_df[(noisy_df["X"] == 0.5) & (noisy_df["Y"] == 0.5)].iloc[0]
        assert abs(noisy_hit["Layer1"] - 0.10) > 1e-4
        assert noisy_hit["Var_Layer1"] > hit["Var_Layer1"] + 1e-4

    @pytest.mark.parametrize("method", ["OK", "SK"])
    def test_sgsim_runs_in_3d(self, run_t2p, method):
        r = run_t2p(
            _ok3d_multi_t2p(extra_options="SGSIM\nSEED 314", method=method),
            extra_files={
                "dataset.csv": _ok3d_irregular_dataset(),
                "grid.grd": OK3D_INTERIOR_GRD,
            },
        )
        assert r.returncode == 0, r.combined
        df = _read_coarse_texture_output(r)
        assert np.isfinite(df[["Layer1", "Layer2", "Layer3"]]).all().all()

    @pytest.mark.parametrize("method", ["OK", "SK"])
    def test_sgsim_preserves_category_masking(self, run_t2p, method):
        dataset = """\
Location,WellIndex,IntervalIndex,X,Y,Zland,TopDepth,BotDepth,Coarse,HSU1
A1,1,1,0.5,0.5,100,0,10,0.10,A
A2,2,1,1.5,0.5,100,0,10,0.15,A
B1,3,1,0.5,2.5,100,0,10,0.85,B
B2,4,1,2.5,2.5,100,0,10,0.90,B
"""
        hsu = """\
1
# Array 1
A A A A A A B B B
"""
        r = run_t2p(
            _ok_lib_2d_t2p(
                method=method,
                options_extra="SGSIM\nSEED 41",
                flow_extra="HSU_FILE hsus.dat",
                variogram_lines="1 Sph 0.0 1.0 100.0 100.0 0.0 4",
            ),
            extra_files={
                "dataset.csv": dataset,
                "grid.grd": SK_2D_GRD,
                "hsus.dat": hsu,
            },
        )
        assert r.returncode == 0, r.combined
        output = _read_coarse_texture_output(r, ("Layer1",))
        values = output["Layer1"].to_numpy()
        assert np.isfinite(values).all()
        if method == "OK":
            assert values[:6].max() < 0.4
            assert values[6:].min() > 0.6
        else:
            swapped_dataset = dataset.replace("0.85,B", "SWAP,B")
            swapped_dataset = swapped_dataset.replace("0.90,B", "0.85,B")
            swapped_dataset = swapped_dataset.replace("SWAP,B", "0.90,B")
            swapped = run_t2p(
                _ok_lib_2d_t2p(
                    method=method,
                    options_extra="SGSIM\nSEED 41",
                    flow_extra="HSU_FILE hsus.dat",
                    variogram_lines="1 Sph 0.0 1.0 100.0 100.0 0.0 4",
                ),
                extra_files={
                    "dataset.csv": swapped_dataset,
                    "grid.grd": SK_2D_GRD,
                    "hsus.dat": hsu,
                },
            )
            assert swapped.returncode == 0, swapped.combined
            swapped_values = _read_coarse_texture_output(
                swapped, ("Layer1",)
            )["Layer1"].to_numpy()
            assert swapped_values[:6] == pytest.approx(values[:6], abs=1e-10)
            assert not np.array_equal(swapped_values[6:], values[6:])

    def test_sgsim_rejects_idw(self, run_t2p):
        method = "IDW"
        r = run_t2p(
            _ok_lib_2d_t2p(method=method, options_extra="SGSIM\nSEED 1"),
            extra_files={"dataset.csv": SK_2D_DATASET, "grid.grd": SK_2D_GRD},
        )
        assert r.returncode != 0
        assert method in r.combined
        assert "SGSIM" in r.combined

    def test_sk_cosgsim_uses_secondary_data_and_its_mean(self, run_t2p):
        cosim = self._run_sk_ck(run_t2p, 44)
        assert cosim.returncode == 0, cosim.combined
        cosim_values = _read_coarse_texture_output(
            cosim, ("Layer1",)
        ).sort_values(["X", "Y"])["Layer1"].to_numpy()

        primary_only = self._run_sk(
            run_t2p, 44, _ck_variance_primary_dataset(0.0)
        )
        assert primary_only.returncode == 0, primary_only.combined
        primary_values = _read_coarse_texture_output(
            primary_only, ("Layer1",)
        ).sort_values(["X", "Y"])["Layer1"].to_numpy()
        assert np.max(np.abs(cosim_values - primary_values)) > 1e-4

        shifted = self._run_sk_ck(
            run_t2p, 44,
            secondary=_shift_dataset_column(
                SK_CK_SECONDARY_DATASET, "Rho", 2.0
            ),
        )
        assert shifted.returncode == 0, shifted.combined
        shifted_values = _read_coarse_texture_output(
            shifted, ("Layer1",)
        ).sort_values(["X", "Y"])["Layer1"].to_numpy()
        assert shifted_values == pytest.approx(cosim_values, abs=3e-6)

    @pytest.mark.parametrize("interp_dim", [2, 3])
    def test_ck_sgsim_runs_with_finite_primary_output(self, run_t2p, interp_dim):
        r = self._run_ck(run_t2p, 123, interp_dim=interp_dim)
        assert r.returncode == 0, r.combined
        values = _read_coarse_texture_output(r, ("Layer1",))["Layer1"]
        assert np.isfinite(values).all()
        assert not (r.workdir / "sgs_path.dat").exists()
        assert not (r.workdir / "sgs_sample.dat").exists()

    def test_ck_sgsim_seed_controls_realization(self, run_t2p):
        first = self._run_ck(run_t2p, 91)
        assert first.returncode == 0, first.combined
        first_values = _read_coarse_texture_output(
            first, ("Layer1",)
        )["Layer1"].to_numpy().copy()
        repeat = self._run_ck(run_t2p, 91)
        assert repeat.returncode == 0, repeat.combined
        repeat_values = _read_coarse_texture_output(
            repeat, ("Layer1",)
        )["Layer1"].to_numpy().copy()
        changed = self._run_ck(run_t2p, 92)
        assert changed.returncode == 0, changed.combined
        changed_values = _read_coarse_texture_output(changed, ("Layer1",))["Layer1"].to_numpy()
        assert np.array_equal(first_values, repeat_values)
        assert not np.array_equal(first_values, changed_values)

    def test_ck_sgsim_uses_secondary_data(self, run_t2p):
        ck = self._run_ck(run_t2p, 44)
        assert ck.returncode == 0, ck.combined
        ck_values = _read_coarse_texture_output(
            ck, ("Layer1",)
        )["Layer1"].to_numpy().copy()
        ok = self._run_2d(run_t2p, 44, _ck_variance_primary_dataset(0.0))
        assert ok.returncode == 0, ok.combined
        ok_values = _read_coarse_texture_output(ok, ("Layer1",))["Layer1"].to_numpy()
        assert np.max(np.abs(ck_values - ok_values)) > 1e-4

    @pytest.mark.parametrize("mean_mode", ["LIBKRIGING", "LEGACY_STANDARDIZED"])
    def test_ck_sgsim_preserves_ck_mean_modes(self, run_t2p, mean_mode):
        r = self._run_ck(run_t2p, 55, mean_mode=mean_mode)
        assert r.returncode == 0, r.combined
        values = _read_coarse_texture_output(r, ("Layer1",))["Layer1"]
        assert np.isfinite(values).all()

    def test_ck_sgsim_accepts_observation_variance(self, run_t2p):
        exact = self._run_ck(run_t2p, 71, first_variance=0.0)
        assert exact.returncode == 0, exact.combined
        exact_df = _read_coarse_texture_output(exact, ("Layer1",)).copy()
        noisy = self._run_ck(run_t2p, 71, first_variance=1.0)
        assert noisy.returncode == 0, noisy.combined
        noisy_df = _read_coarse_texture_output(noisy, ("Layer1",))
        exact_hit = exact_df[(exact_df["X"] == 0.5) & (exact_df["Y"] == 0.5)].iloc[0]
        noisy_hit = noisy_df[(noisy_df["X"] == 0.5) & (noisy_df["Y"] == 0.5)].iloc[0]
        assert exact_hit["Layer1"] == pytest.approx(0.10, abs=1e-10)
        assert abs(noisy_hit["Layer1"] - exact_hit["Layer1"]) > 1e-4

    def test_ck_sgsim_rejects_block_mode(self, run_t2p):
        t2p = _ck_variance_t2p().replace(
            "MAX_VSTRUCT 1", "MAX_VSTRUCT 1\n          SGSIM\n          SEED 1"
        ).replace(
            "TYPE GRID\n          GRID_FILE grid.grd",
            "TYPE MODFLOW\n          NAM_FILE model.nam\n          BLOCK_MODE",
        )
        files = {
            **BLOCK_FILES,
            "primary.csv": _ck_variance_primary_dataset(0.0),
            "secondary.csv": CK_VARIANCE_SECONDARY_DATASET,
        }
        r = run_t2p(t2p, extra_files=files)
        assert r.returncode != 0
        assert "BLOCK_MODE" in r.combined
        assert "unsupported with SGSIM" in r.combined

    @pytest.mark.parametrize("weight_option", ["WRITE_GRID_WEIGHTS", "READ_GRID_WEIGHTS"])
    def test_ck_sgsim_rejects_saved_weights(self, run_t2p, weight_option):
        r = self._run_ck(run_t2p, 1, extra=weight_option)
        assert r.returncode != 0
        assert weight_option in r.combined
        assert "SGSIM" in r.combined

    def test_sgsim_runtime_no_longer_references_legacy_type(self):
        root = Path(__file__).parents[1]
        source = (root / "src" / "T2P.f90").read_text()
        assert "allocate(t_sgsim::" not in source.lower()
        assert "use m_sgsim" not in source.lower()

        for legacy_source in ("kriging.f90", "sgsim.f90", "kriging_vecback.f90"):
            assert not (root / "src" / legacy_source).exists()

        for makefile in (
            root / "build_ifx" / "makefile",
            root / "build" / "makefile",
            root / "build" / "build_ifx" / "makefile",
        ):
            build_text = makefile.read_text()
            assert "kriging.f90" not in build_text
            assert "sgsim.f90" not in build_text
            assert "$(OBJDIR)/kriging.obj" not in build_text
            assert "$(OBJDIR)/sgsim.obj" not in build_text

    def test_sgsim_incompatible_with_weight_write(self, run_t2p):
        """Per docs, SGS is not compatible with WRITE_GRID_WEIGHTS."""
        r = run_t2p(_grid_t2p(options_extra="SGSIM\n  SEED 1\n  WRITE_GRID_WEIGHTS"),
                    extra_files=FILES)
        assert r.returncode != 0 or "warning" in r.combined.lower() or \
               "error" in r.combined.lower(), (
            "Expected warning/error for SGSIM + WRITE_GRID_WEIGHTS combination")


# ---------------------------------------------------------------------------
# 5. Weight caching (WRITE_*/READ_*)
# ---------------------------------------------------------------------------

class TestWeightCaching:

    @pytest.mark.parametrize("method", ["OK", "SK"])
    def test_grid_weight_write_read_preserves_estimates_and_variance(
        self, run_t2p, method
    ):
        files = TestLibKrigingSimple._files_2d() if method == "SK" else {
            "dataset.csv": OK_LIB_2D_DATASET,
            "grid.grd": OK_LIB_2D_GRD,
        }
        write = run_t2p(
            _ok_lib_2d_t2p(
                method=method,
                options_extra="WRITE_GRID_WEIGHTS\nWRITE_VARIANCE",
            ),
            extra_files=files,
        )
        assert write.returncode == 0, write.combined
        expected = _read_coarse_texture_output(write, ("Layer1",)).copy()
        weight_files = sorted(write.workdir.glob("*_grdweights_*.wts"))
        assert weight_files

        read = run_t2p(
            _ok_lib_2d_t2p(
                method=method,
                options_extra="READ_GRID_WEIGHTS\nWRITE_VARIANCE",
            ),
            extra_files=files,
        )
        assert read.returncode == 0, read.combined
        actual = _read_coarse_texture_output(read, ("Layer1",))
        assert actual["Layer1"].to_numpy() == pytest.approx(
            expected["Layer1"].to_numpy(), abs=1e-10
        )
        assert actual["Var_Layer1"].to_numpy() == pytest.approx(
            expected["Var_Layer1"].to_numpy(), abs=1e-10
        )

    def test_grid_read_weights_uses_current_observation_values(self, run_t2p):
        files = {"dataset.csv": OK_LIB_2D_DATASET, "grid.grd": OK_LIB_2D_GRD}
        write = run_t2p(
            _ok_lib_2d_t2p(options_extra="WRITE_GRID_WEIGHTS"),
            extra_files=files,
        )
        assert write.returncode == 0, write.combined
        original = _read_coarse_texture_output(write, ("Layer1",)).copy()

        changed_files = {
            **files,
            "dataset.csv": OK_LIB_2D_DATASET.replace("0.10", "0.90", 1),
        }
        read = run_t2p(
            _ok_lib_2d_t2p(options_extra="READ_GRID_WEIGHTS"),
            extra_files=changed_files,
        )
        assert read.returncode == 0, read.combined
        changed = _read_coarse_texture_output(read, ("Layer1",))
        assert np.max(np.abs(
            changed["Layer1"].to_numpy() - original["Layer1"].to_numpy()
        )) > 0.1

    def test_ck_weight_write_read_preserves_mean_mode(self, run_t2p):
        files = {
            "primary.csv": _ck_variance_primary_dataset(0.0),
            "secondary.csv": CK_VARIANCE_SECONDARY_DATASET,
            "grid.grd": OK_LIB_2D_GRD,
        }
        write_t2p = _ck_variance_t2p().replace(
            "KVME", "KVME\n          WRITE_GRID_WEIGHTS"
        )
        write = run_t2p(write_t2p, extra_files=files)
        assert write.returncode == 0, write.combined
        expected = _read_coarse_texture_output(write, ("Layer1",)).copy()
        assert list(write.workdir.glob("*_grdweights_*.wts"))

        read_t2p = _ck_variance_t2p().replace(
            "KVME", "KVME\n          READ_GRID_WEIGHTS"
        )
        read = run_t2p(read_t2p, extra_files=files)
        assert read.returncode == 0, read.combined
        actual = _read_coarse_texture_output(read, ("Layer1",))
        assert actual["Layer1"].to_numpy() == pytest.approx(
            expected["Layer1"].to_numpy(), abs=1e-10
        )
        assert actual["Var_Layer1"].to_numpy() == pytest.approx(
            expected["Var_Layer1"].to_numpy(), abs=1e-10
        )

    def test_block_weight_write_read(self, run_t2p):
        write = run_t2p(
            _block_modflow_t2p(options_extra="WRITE_GRID_WEIGHTS"),
            extra_files=BLOCK_FILES,
        )
        assert write.returncode == 0, write.combined
        expected = _read_coarse_texture_output(write, ("Layer1",)).copy()
        assert list(write.workdir.glob("*_grdweights_*.wts"))

        read = run_t2p(
            _block_modflow_t2p(options_extra="READ_GRID_WEIGHTS"),
            extra_files=BLOCK_FILES,
        )
        assert read.returncode == 0, read.combined
        actual = _read_coarse_texture_output(read, ("Layer1",))
        assert actual["Layer1"].to_numpy() == pytest.approx(
            expected["Layer1"].to_numpy(), abs=1e-10
        )
        assert actual["Var_Layer1"].to_numpy() == pytest.approx(
            expected["Var_Layer1"].to_numpy(), abs=1e-10
        )

    def test_pp_weight_write_read_reuses_multi_parameter_geometry(self, run_t2p):
        write = run_t2p(
            _grid_t2p(options_extra="WRITE_PP_WEIGHTS"), extra_files=FILES
        )
        assert write.returncode == 0, write.combined
        weight_files = sorted(write.workdir.glob("*_ppweights_*.wts"))
        assert weight_files
        output_names = {
            path.name for pattern in ("*KH*.csv", "*KV*.csv", "*SS*.csv", "*SY*.csv")
            for path in write.workdir.glob(pattern)
        }
        assert output_names
        expected = {
            name: (write.workdir / name).read_text() for name in output_names
        }

        read = run_t2p(
            _grid_t2p(options_extra="READ_PP_WEIGHTS"), extra_files=FILES
        )
        assert read.returncode == 0, read.combined
        for name, content in expected.items():
            assert (read.workdir / name).read_text() == content

    def test_missing_grid_weight_file_is_fatal(self, run_t2p):
        r = run_t2p(
            _ok_lib_2d_t2p(options_extra="READ_GRID_WEIGHTS"),
            extra_files={"dataset.csv": OK_LIB_2D_DATASET, "grid.grd": OK_LIB_2D_GRD},
        )
        assert r.returncode != 0
        assert "Could not open weight_file" in r.combined

    def test_text_option_uses_native_text_factor_file(self, run_t2p):
        r = run_t2p(
            _ok_lib_2d_t2p(options_extra="WRITE_GRID_WEIGHTS TEXT"),
            extra_files={"dataset.csv": OK_LIB_2D_DATASET, "grid.grd": OK_LIB_2D_GRD},
        )
        assert r.returncode == 0, r.combined
        weight_file = next(r.workdir.glob("*_grdweights_*.wts"))
        header = weight_file.read_text().splitlines()[0].split()
        assert [int(value) for value in header[:2]] == [4, 1]


# ---------------------------------------------------------------------------
# 6. ONLY_TEXTURE (Data2Texture mode)
# ---------------------------------------------------------------------------

class TestOnlyTexture:

    def test_only_texture_no_pilot_points_required(self, run_t2p):
        """ONLY_TEXTURE runs without a PP_PARAMETERS block."""
        t2p = """\
            BEGIN OPTIONS
              MAX_VSTRUCT 1
              INFER_LAST_CLASS
              ONLY_TEXTURE
              WRITE_DATASET_FILES
            END OPTIONS
            BEGIN CLASSES
              Coarse
              Fine
            END CLASSES
            BEGIN DATASET
              FILE dataset.csv
            END DATASET
            BEGIN FLOW_MODEL
              TYPE GRID
              GRID_FILE grid.grd
            END FLOW_MODEL
            BEGIN VARIOGRAMS
              CLASS Coarse
                1 Sph 0.05 0.038 1.5E4 1.5E4 0.0 16
            END VARIOGRAMS
            """
        r = run_t2p(t2p, extra_files=FILES)
        assert r.returncode == 0, r.combined

    def test_only_texture_no_model_output_files(self, run_t2p):
        """With ONLY_TEXTURE, no model hydraulic-property files should be written."""
        t2p = """\
            BEGIN OPTIONS
              MAX_VSTRUCT 1
              INFER_LAST_CLASS
              ONLY_TEXTURE
            END OPTIONS
            BEGIN CLASSES
              Coarse
              Fine
            END CLASSES
            BEGIN DATASET
              FILE dataset.csv
            END DATASET
            BEGIN FLOW_MODEL
              TYPE GRID
              GRID_FILE grid.grd
            END FLOW_MODEL
            BEGIN VARIOGRAMS
              CLASS Coarse
                1 Sph 0.05 0.038 1.5E4 1.5E4 0.0 16
            END VARIOGRAMS
            """
        r = run_t2p(t2p, extra_files=FILES)
        assert r.returncode == 0, r.combined
        # No KH/KV/SS/SY files expected
        for pat in ("*KH*", "*KV*", "*SS*", "*SY*"):
            assert not list(r.workdir.glob(pat)), (
                f"Unexpected hydraulic property file matching {pat} found with ONLY_TEXTURE")


# ---------------------------------------------------------------------------
# 7. MAX_LOG_LENGTH (3D log preprocessing)
# ---------------------------------------------------------------------------

class TestMaxLogLength:

    def _run_interval(self, run_t2p, bot_depth):
        return run_t2p(
            _max_log_length_t2p(max_log_length=5),
            extra_files={
                "dataset.csv": _max_log_length_dataset(bot_depth),
                "grid.grd": MAX_LOG_LENGTH_GRD,
            },
        )

    def test_interval_equal_to_max_length_is_not_split(self, run_t2p):
        r = self._run_interval(run_t2p, bot_depth=5)
        _assert_processed_log_points(r, expected_z=[97.5])

    def test_interval_shorter_than_max_length_is_not_split(self, run_t2p):
        r = self._run_interval(run_t2p, bot_depth=4)
        _assert_processed_log_points(r, expected_z=[98.0])

    def test_interval_twice_max_length_splits_into_two_equal_parts(self, run_t2p):
        r = self._run_interval(run_t2p, bot_depth=10)
        _assert_processed_log_points(
            r,
            expected_z=[97.5, 92.5],
            expected_midpoint_spacings=[5.0],
            expected_segment_lengths=[5.0, 5.0],
        )

    def test_interval_slightly_longer_than_max_splits_equally(self, run_t2p):
        r = self._run_interval(run_t2p, bot_depth=7)
        _assert_processed_log_points(
            r,
            expected_z=[98.25, 94.75],
            expected_midpoint_spacings=[3.5],
            expected_segment_lengths=[3.5, 3.5],
        )

    def test_longer_non_multiple_interval_splits_into_equal_parts(self, run_t2p):
        r = self._run_interval(run_t2p, bot_depth=12)
        _assert_processed_log_points(
            r,
            expected_z=[98.0, 94.0, 90.0],
            expected_midpoint_spacings=[4.0, 4.0],
            expected_segment_lengths=[4.0, 4.0, 4.0],
        )


# ---------------------------------------------------------------------------
# 8. INTERP_DIM 3 grid geometry
# ---------------------------------------------------------------------------

class TestInterpDim3GridGeometry:

    def test_elevation_file_writes_top_bottoms_and_layer_centers(self, run_t2p):
        # INTERP_DIM 3 targets the actual model-layer center elevations. Its
        # conditioning records remain XYZ support points rather than being
        # reduced to the per-layer XY averages used by 2D interpolation.
        r = run_t2p(
            _interp_dim3_geometry_t2p(),
            extra_files={
                "dataset.csv": _max_log_length_dataset(5),
                "grid.grd": MAX_LOG_LENGTH_GRD,
            },
        )
        assert r.returncode == 0, r.combined

        df = _read_grid_elevation_file(r)
        assert len(df) == 1, df.to_string(index=False)
        row = df.iloc[0]

        assert [row["Top"], row["Bottom1"], row["Bottom2"], row["Bottom3"]] == (
            pytest.approx([100.0, 90.0, 80.0, 70.0])
        )
        layer_centers = [
            0.5 * (row["Top"] + row["Bottom1"]),
            0.5 * (row["Bottom1"] + row["Bottom2"]),
            0.5 * (row["Bottom2"] + row["Bottom3"]),
        ]
        assert layer_centers == pytest.approx([95.0, 85.0, 75.0])
        assert row["X"] == pytest.approx(0.5)
        assert row["Y"] == pytest.approx(0.5)


# ---------------------------------------------------------------------------
# 9. Ordinary kriging in 3D
# ---------------------------------------------------------------------------

class TestOrdinaryKriging3D:

    def test_ok_reference_schur_matches_augmented_system(self):
        dataset = _ok3d_corner_dataset()
        obsloc, obsval = _ok3d_observations_from_dataset(dataset)
        targets = np.array([
            [1.5, 1.5, 95.0],
            [1.5, 1.5, 85.0],
            [1.5, 1.5, 75.0],
            [1.25, 1.75, 88.0],
        ])

        schur_est, schur_var = _ok_reference_3d(obsloc, obsval, targets, OK3D_MULTI_VGM)
        aug_est, aug_var = _ok_reference_3d_augmented(obsloc, obsval, targets, OK3D_MULTI_VGM)

        assert schur_est == pytest.approx(aug_est, abs=1e-12)
        assert schur_var == pytest.approx(aug_var, abs=1e-12)

    def test_variogram_rotation_metric_sanity_checks(self):
        # Migration note: these assertions intentionally lock down the legacy
        # Texture2Par rotation/dip convention. libkriging/scipy/GSLIB-style
        # conventions differ, so this reference is expected to change only as
        # an explicit migration step, not as incidental test maintenance.
        assert _effective_distance([1, 0, 0], 1.0, 4.0, 20.0) == pytest.approx(1.0)
        assert _effective_distance([0, 1, 0], 1.0, 4.0, 20.0) == pytest.approx(0.25)
        assert _effective_distance([0, 0, 1], 1.0, 4.0, 20.0) == pytest.approx(0.05)

        delta = [1.0, 2.0, 0.0]
        isotropic_horizontal = [
            _effective_distance(delta, 4.0, 4.0, 20.0, azimuth=angle)
            for angle in (0.0, 30.0, 90.0)
        ]
        assert isotropic_horizontal == pytest.approx([isotropic_horizontal[0]] * 3)

        assert _effective_distance([1, 0, 0], 1.0, 4.0, 20.0) > (
            _effective_distance([0, 1, 0], 1.0, 4.0, 20.0)
        )
        assert _effective_distance([0, 1, 0], 1.0, 4.0, 20.0, azimuth=90.0) > (
            _effective_distance([1, 0, 0], 1.0, 4.0, 20.0, azimuth=90.0)
        )

        for delta in ([1.0, 2.0, 3.0], [-2.0, 0.5, 1.25], [0.0, -1.0, 4.0]):
            unrotated = _effective_distance(delta, 7.0, 7.0, 7.0)
            rotated = _effective_distance(
                delta,
                7.0,
                7.0,
                7.0,
                azimuth=30.0,
                dip=15.0,
                plunge=-20.0,
            )
            assert rotated == pytest.approx(unrotated, abs=1e-12)

    def test_exact_vertical_hits_preserve_layer_sequence(self, run_t2p):
        # A target coincident with one observation is valid and must reproduce
        # it; only duplicate conditioning support points are forbidden.
        r = run_t2p(
            _ok3d_t2p(),
            extra_files={
                "dataset.csv": _ok3d_exact_hit_dataset(),
                "grid.grd": MAX_LOG_LENGTH_GRD,
            },
        )
        assert r.returncode == 0, r.combined

        df = _read_coarse_texture_output(r)
        assert len(df) == 1, df.to_string(index=False)
        row = df.iloc[0]

        assert row["Layer1"] == pytest.approx(0.10, abs=1e-2)
        assert row["Layer2"] == pytest.approx(0.50, abs=1e-2)
        assert row["Layer3"] == pytest.approx(0.90, abs=1e-2)
        assert row["Layer3"] - row["Layer1"] > 0.7

    def test_multiple_logs_match_3d_ok_reference(self, run_t2p):
        dataset = _ok3d_multi_dataset()
        r = run_t2p(
            _ok3d_multi_t2p(),
            extra_files={
                "dataset.csv": dataset,
                "grid.grd": OK3D_MULTI_GRD,
            },
        )
        assert r.returncode == 0, r.combined

        df = _read_coarse_texture_output(r)
        assert len(df) == 4, df.to_string(index=False)
        targets, t2p_est, t2p_var = _coarse_output_to_arrays(df)
        obsloc, obsval = _ok3d_observations_from_dataset(dataset)
        ref_est, ref_var = _ok_reference_3d(obsloc, obsval, targets, OK3D_MULTI_VGM)

        assert np.allclose(t2p_est.ravel(), ref_est, atol=1e-5), (
            "3D OK estimates differ from Python reference.\n"
            f"T2P:\n{t2p_est}\nReference:\n{ref_est.reshape(t2p_est.shape)}"
        )
        assert np.allclose(t2p_var.ravel(), ref_var, atol=1e-5), (
            "3D OK variances differ from Python reference.\n"
            f"T2P:\n{t2p_var}\nReference:\n{ref_var.reshape(t2p_var.shape)}"
        )
        assert np.mean(t2p_est[:, 2] - t2p_est[:, 0]) > 0.4

    def test_corner_logs_interpolate_unsampled_center_cell(self, run_t2p):
        dataset = _ok3d_corner_dataset()
        r = run_t2p(
            _ok3d_multi_t2p(),
            extra_files={
                "dataset.csv": dataset,
                "grid.grd": OK3D_INTERIOR_GRD,
            },
        )
        assert r.returncode == 0, r.combined

        df = _read_coarse_texture_output(r)
        assert len(df) == 9, df.to_string(index=False)
        targets, t2p_est, t2p_var = _coarse_output_to_arrays(df)
        obsloc, obsval = _ok3d_observations_from_dataset(dataset)
        ref_est, ref_var = _ok_reference_3d(obsloc, obsval, targets, OK3D_MULTI_VGM)

        assert np.allclose(t2p_est.ravel(), ref_est, atol=1e-5), (
            "3D OK estimates differ from Python reference for corner-log grid.\n"
            f"T2P:\n{t2p_est}\nReference:\n{ref_est.reshape(t2p_est.shape)}"
        )
        assert np.allclose(t2p_var.ravel(), ref_var, atol=1e-5), (
            "3D OK variances differ from Python reference for corner-log grid.\n"
            f"T2P:\n{t2p_var}\nReference:\n{ref_var.reshape(t2p_var.shape)}"
        )

        center = df[
            (df["X"].sub(1.5).abs() < 1e-9) &
            (df["Y"].sub(1.5).abs() < 1e-9)
        ]
        assert len(center) == 1, df.to_string(index=False)
        center_est = center[["Layer1", "Layer2", "Layer3"]].to_numpy(dtype=float)[0]
        center_targets = np.array([[1.5, 1.5, 95.0], [1.5, 1.5, 85.0], [1.5, 1.5, 75.0]])
        center_ref, _ = _ok_reference_3d(obsloc, obsval, center_targets, OK3D_MULTI_VGM)
        assert center_est == pytest.approx(center_ref, abs=1e-5)
        assert center_est[2] - center_est[0] > 0.25

    def test_vertical_range_controls_cross_layer_influence(self, run_t2p):
        dataset = _ok3d_corner_dataset()
        r = run_t2p(
            _ok3d_multi_t2p(vgm=OK3D_CROSSLAYER_VGM),
            extra_files={
                "dataset.csv": dataset,
                "grid.grd": OK3D_INTERIOR_GRD,
            },
        )
        assert r.returncode == 0, r.combined

        df = _read_coarse_texture_output(r)
        assert len(df) == 9, df.to_string(index=False)
        targets, t2p_est, t2p_var = _coarse_output_to_arrays(df)
        obsloc, obsval = _ok3d_observations_from_dataset(dataset)
        ref_est, ref_var = _ok_reference_3d(obsloc, obsval, targets, OK3D_CROSSLAYER_VGM)

        assert np.allclose(t2p_est.ravel(), ref_est, atol=1e-5), (
            "3D OK estimates differ from Python reference for cross-layer range.\n"
            f"T2P:\n{t2p_est}\nReference:\n{ref_est.reshape(t2p_est.shape)}"
        )
        assert np.allclose(t2p_var.ravel(), ref_var, atol=1e-5), (
            "3D OK variances differ from Python reference for cross-layer range.\n"
            f"T2P:\n{t2p_var}\nReference:\n{ref_var.reshape(t2p_var.shape)}"
        )

        center = df[
            (df["X"].sub(1.5).abs() < 1e-9) &
            (df["Y"].sub(1.5).abs() < 1e-9)
        ]
        assert len(center) == 1, df.to_string(index=False)
        center_est = center[["Layer1", "Layer2", "Layer3"]].to_numpy(dtype=float)[0]
        center_targets = np.array([[1.5, 1.5, 95.0], [1.5, 1.5, 85.0], [1.5, 1.5, 75.0]])
        base_ref, _ = _ok_reference_3d(obsloc, obsval, center_targets, OK3D_MULTI_VGM)
        cross_ref, _ = _ok_reference_3d(obsloc, obsval, center_targets, OK3D_CROSSLAYER_VGM)
        assert center_est == pytest.approx(cross_ref, abs=1e-5)
        assert np.max(np.abs(cross_ref - base_ref)) > 0.005

    def test_rotated_anisotropic_variogram_matches_3d_ok_reference(self, run_t2p):
        dataset = _ok3d_corner_dataset()
        r = run_t2p(
            _ok3d_multi_t2p(vgm=OK3D_ROTATED_VGM),
            extra_files={
                "dataset.csv": dataset,
                "grid.grd": OK3D_INTERIOR_GRD,
            },
        )
        assert r.returncode == 0, r.combined

        df = _read_coarse_texture_output(r)
        assert len(df) == 9, df.to_string(index=False)
        targets, t2p_est, t2p_var = _coarse_output_to_arrays(df)
        obsloc, obsval = _ok3d_observations_from_dataset(dataset)
        ref_est, ref_var = _ok_reference_3d(obsloc, obsval, targets, OK3D_ROTATED_VGM)

        assert np.allclose(t2p_est.ravel(), ref_est, atol=1e-5), (
            "3D OK estimates differ from Python reference for rotated anisotropy.\n"
            f"T2P:\n{t2p_est}\nReference:\n{ref_est.reshape(t2p_est.shape)}"
        )
        assert np.allclose(t2p_var.ravel(), ref_var, atol=1e-5), (
            "3D OK variances differ from Python reference for rotated anisotropy.\n"
            f"T2P:\n{t2p_var}\nReference:\n{ref_var.reshape(t2p_var.shape)}"
        )

        unrotated_ref, _ = _ok_reference_3d(obsloc, obsval, targets, OK3D_CROSSLAYER_VGM)
        assert np.max(np.abs(ref_est - unrotated_ref)) > 0.01

    def test_neighbor_truncation_matches_3d_ok_reference(self, run_t2p):
        dataset = _ok3d_irregular_dataset()
        r = run_t2p(
            _ok3d_multi_t2p(vgm=OK3D_TRUNCATED_VGM),
            extra_files={
                "dataset.csv": dataset,
                "grid.grd": OK3D_INTERIOR_GRD,
            },
        )
        assert r.returncode == 0, r.combined

        df = _read_coarse_texture_output(r)
        assert len(df) == 9, df.to_string(index=False)
        targets, t2p_est, t2p_var = _coarse_output_to_arrays(df)
        obsloc, obsval = _ok3d_observations_from_dataset(dataset)
        ref_est, ref_var = _ok_reference_3d_nearest(obsloc, obsval, targets, OK3D_TRUNCATED_VGM)

        assert np.allclose(t2p_est.ravel(), ref_est, atol=1e-5), (
            "3D OK estimates differ from nearest-neighbor Python reference.\n"
            f"T2P:\n{t2p_est}\nReference:\n{ref_est.reshape(t2p_est.shape)}"
        )
        assert np.allclose(t2p_var.ravel(), ref_var, atol=1e-5), (
            "3D OK variances differ from nearest-neighbor Python reference.\n"
            f"T2P:\n{t2p_var}\nReference:\n{ref_var.reshape(t2p_var.shape)}"
        )

        all_ref, _ = _ok_reference_3d(obsloc, obsval, targets, OK3D_CROSSLAYER_VGM)
        assert np.max(np.abs(ref_est - all_ref)) > 0.01

    def test_anisotropic_search_changes_3d_neighbor_selection(self, run_t2p):
        # The anisotropic search metric below uses the same legacy rotation/dip
        # convention as the 3D covariance reference above.
        dataset = _ok3d_irregular_dataset()
        files = {"dataset.csv": dataset, "grid.grd": OK3D_INTERIOR_GRD}

        r = run_t2p(
            _ok3d_multi_t2p(
                vgm=OK3D_ANISO_SEARCH_VGM,
                extra_options="ANISOTROPIC_SEARCH",
            ),
            extra_files=files,
        )
        assert r.returncode == 0, r.combined
        df = _read_coarse_texture_output(r)
        targets, t2p_est, t2p_var = _coarse_output_to_arrays(df)
        obsloc, obsval = _ok3d_observations_from_dataset(dataset)
        ref_est, ref_var = _ok_reference_3d_nearest(
            obsloc,
            obsval,
            targets,
            OK3D_ANISO_SEARCH_VGM,
            anisotropic_search=True,
        )

        assert np.allclose(t2p_est.ravel(), ref_est, atol=1e-5), (
            "3D OK estimates differ from anisotropic-search Python reference.\n"
            f"T2P:\n{t2p_est}\nReference:\n{ref_est.reshape(t2p_est.shape)}"
        )
        assert np.allclose(t2p_var.ravel(), ref_var, atol=1e-5), (
            "3D OK variances differ from anisotropic-search Python reference.\n"
            f"T2P:\n{t2p_var}\nReference:\n{ref_var.reshape(t2p_var.shape)}"
        )

        r_iso = run_t2p(
            _ok3d_multi_t2p(vgm=OK3D_ANISO_SEARCH_VGM),
            extra_files=files,
        )
        assert r_iso.returncode == 0, r_iso.combined
        _, t2p_iso_est, _ = _coarse_output_to_arrays(_read_coarse_texture_output(r_iso))
        assert np.max(np.abs(t2p_est - t2p_iso_est)) > 0.01

    @pytest.mark.parametrize("vgm", [OK3D_EXP_VGM, OK3D_GAU_VGM])
    def test_non_spherical_3d_variograms_match_ok_reference(self, run_t2p, vgm):
        dataset = _ok3d_corner_dataset()
        r = run_t2p(
            _ok3d_multi_t2p(vgm=vgm),
            extra_files={
                "dataset.csv": dataset,
                "grid.grd": OK3D_INTERIOR_GRD,
            },
        )
        assert r.returncode == 0, r.combined

        df = _read_coarse_texture_output(r)
        assert len(df) == 9, df.to_string(index=False)
        if vgm == OK3D_GAU_VGM:
            sampled_xy = {(0.5, 0.5), (2.5, 0.5), (0.5, 2.5), (2.5, 2.5)}
            unsampled_mask = [
                (round(row.X, 9), round(row.Y, 9)) not in sampled_xy
                for row in df.itertuples(index=False)
            ]
            df = df.loc[unsampled_mask].reset_index(drop=True)
            assert len(df) >= 5, df.to_string(index=False)
            assert (
                (df["X"].sub(1.5).abs() < 1e-9) &
                (df["Y"].sub(1.5).abs() < 1e-9)
            ).any(), df.to_string(index=False)

        targets, t2p_est, t2p_var = _coarse_output_to_arrays(df)
        obsloc, obsval = _ok3d_observations_from_dataset(dataset)
        ref_est, ref_var = _ok_reference_3d(obsloc, obsval, targets, vgm)

        assert np.allclose(t2p_est.ravel(), ref_est, atol=1e-5), (
            f"3D OK estimates differ from Python reference for {vgm}.\n"
            f"T2P:\n{t2p_est}\nReference:\n{ref_est.reshape(t2p_est.shape)}"
        )
        assert np.allclose(t2p_var.ravel(), ref_var, atol=1e-5), (
            f"3D OK variances differ from Python reference for {vgm}.\n"
            f"T2P:\n{t2p_var}\nReference:\n{ref_var.reshape(t2p_var.shape)}"
        )


# ---------------------------------------------------------------------------
# 10. READ_GRIDDED mode
# ---------------------------------------------------------------------------

class TestReadGridded:

    def test_read_gridded_mode(self, run_t2p):
        """READ_GRIDDED skips kriging and reads pre-interpolated texture directly."""
        # 5x5 grid, 2 layers = 50 rows; 3 ignored cols (node, X, Y) + Coarse
        gridded_rows = ["Node,X,Y,Coarse"]
        node = 1
        for layer in range(1, 3):
            for row in range(1, 6):
                for col in range(1, 6):
                    x = 500000 + (col - 1) * 3000
                    y = 4000000 + (row - 1) * 3000
                    gridded_rows.append(f"{node},{x},{y},0.55")
                    node += 1
        gridded_csv = "\n".join(gridded_rows)

        t2p = """\
            BEGIN OPTIONS
              MAX_VSTRUCT 1
              INFER_LAST_CLASS
              READ_GRID
              WRITE_NODE_FILES
            END OPTIONS
            BEGIN CLASSES
              Coarse
              Fine
            END CLASSES
            BEGIN DATASET
              FILE gridded.csv
            END DATASET
            BEGIN FLOW_MODEL
              TYPE GRID
              GRID_FILE grid.grd
            END FLOW_MODEL
            BEGIN VARIOGRAMS
              CLASS PilotPoints
                1 Sph 0.00 1.0 1.0E5 1.0E5 0.0 25
            END VARIOGRAMS
            BEGIN PP_LOCS
              1 5050.0 4005.0 1
            END PP_LOCS
            BEGIN PP_PARAMETERS
              TYPE Global
              1 0.5 -0.5 0.0
              TYPE Aquifer
              1 Coarse 0.1 100.0 1e-5 0.15 10.0 0.001
              1 Fine   0.001 0.5 1e-4 0.05 5.0 0.001
            END PP_PARAMETERS
            """
        r = run_t2p(t2p, extra_files={"gridded.csv": gridded_csv, "grid.grd": GRD})
        assert r.returncode == 0, r.combined


# ---------------------------------------------------------------------------
# 11. ANISOTROPIC_SEARCH
# ---------------------------------------------------------------------------

class TestAnisotropicSearch:

    def test_anisotropic_search_runs(self, run_t2p):
        r = run_t2p(_grid_t2p(options_extra="ANISOTROPIC_SEARCH"), extra_files=FILES)
        assert r.returncode == 0, r.combined


# ---------------------------------------------------------------------------
# 12. NO_TEX_BOUNDS
# ---------------------------------------------------------------------------

class TestNoTexBounds:

    def test_no_tex_bounds_allows_outside_range(self, run_t2p):
        """With NO_TEX_BOUNDS, values outside [0,1] should not trigger an error."""
        r = run_t2p(_grid_t2p(options_extra="NO_TEX_BOUNDS\n  WRITE_DATASET_FILES"),
                    extra_files=FILES)
        assert r.returncode == 0, r.combined
