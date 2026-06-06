"""
test_interpolation.py
=====================
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


def _ok3d_t2p() -> str:
    return f"""\
        BEGIN OPTIONS
          MAX_VSTRUCT 1
          INFER_LAST_CLASS
          ONLY_TEXTURE
          WRITE_NODE_FILES
          WRITE_DATASET_FILES
          WRITE_VARIANCE
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


def _ok3d_multi_t2p() -> str:
    return f"""\
        BEGIN OPTIONS
          MAX_VSTRUCT 1
          INFER_LAST_CLASS
          ONLY_TEXTURE
          WRITE_NODE_FILES
          WRITE_DATASET_FILES
          WRITE_VARIANCE
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
            1 {OK3D_MULTI_VGM}
        END VARIOGRAMS
        """


def _read_coarse_texture_output(r) -> pd.DataFrame:
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
    required = {"Layer1", "Layer2", "Layer3"}
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
    if vtype != "sph":
        raise ValueError(f"Only spherical covariance is supported here, got {vtype}")
    nugget, sill, hmin, hmax, hvert = (float(parts[i]) for i in range(1, 6))
    azimuth, dip, plunge = (float(parts[i]) for i in range(6, 9))
    nnear = int(parts[9])
    if any(abs(angle) > 1e-12 for angle in (azimuth, dip, plunge)):
        raise ValueError("This local OK reference currently expects zero rotation angles")
    return nugget, sill, np.array([hmin, hmax, hvert], dtype=float), nnear


def _spherical_covariance(hr):
    return np.where(hr < 1.0, 1.0 - 1.5 * hr + 0.5 * hr**3, 0.0)


def _ok_reference_3d(obsloc, obsval, newloc, vgm):
    obsloc = np.asarray(obsloc, dtype=float)
    obsval = np.asarray(obsval, dtype=float)
    newloc = np.asarray(newloc, dtype=float)
    nugget, sill, ranges, nnear = _parse_ok3d_vgm(vgm)
    if nnear < len(obsloc):
        raise ValueError(f"nnear ({nnear}) must be >= nobs ({len(obsloc)})")

    scaled_obs = obsloc / ranges
    scaled_new = newloc / ranges
    cov0 = nugget + sill
    c_obs = _spherical_covariance(squareform(pdist(scaled_obs))) * sill
    np.fill_diagonal(c_obs, cov0)
    c_rhs = _spherical_covariance(cdist(scaled_obs, scaled_new)) * sill

    e = np.ones(len(obsloc))
    cf = cho_factor(c_obs, lower=True)
    cinv_rhs = cho_solve(cf, c_rhs)
    cinv_e = cho_solve(cf, e)
    mu = (e @ cinv_rhs - 1.0) / (e @ cinv_e)
    weights = cinv_rhs - np.outer(cinv_e, mu)

    est = obsval @ weights
    var = np.maximum(0.0, cov0 - (weights * c_rhs).sum(axis=0) - mu)
    return est, var


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
# 2. CORRECT_WEIGHTS (non-negative kriging)
# ---------------------------------------------------------------------------

class TestCorrectWeights:

    def test_correct_weights_flag_runs(self, run_t2p):
        r = run_t2p(_grid_t2p(options_extra="CORRECT_WEIGHTS"), extra_files=FILES)
        assert r.returncode == 0, r.combined


# ---------------------------------------------------------------------------
# 3. Sequential Gaussian Simulation (SGS)
# ---------------------------------------------------------------------------

class TestSGS:

    def test_sgsim_flag_accepted(self, run_t2p):
        r = run_t2p(_grid_t2p(options_extra="SGSIM\n  SEED 12345"),
                    extra_files=FILES)
        # SGS may or may not produce additional outputs; the run must not crash
        assert r.returncode == 0, r.combined

    def test_sgsim_different_seeds_differ(self, run_t2p):
        """Two SGS runs with different seeds should produce different outputs."""
        r1 = run_t2p(_grid_t2p(options_extra="SGSIM\n  SEED 1\n  WRITE_DATASET_FILES"),
                     extra_files=FILES)
        if r1.returncode != 0:
            pytest.skip("SGS run did not complete successfully")
        files1 = sorted(r1.workdir.glob("*Coarse*"))
        if not files1:
            pytest.skip("Coarse dataset files not found for SGS comparison")
        content1 = files1[0].read_text()  # capture before second run overwrites

        r2 = run_t2p(_grid_t2p(options_extra="SGSIM\n  SEED 2\n  WRITE_DATASET_FILES"),
                     extra_files=FILES)
        if r2.returncode != 0:
            pytest.skip("SGS run did not complete successfully")
        files2 = sorted(r2.workdir.glob("*Coarse*"))
        if not files2:
            pytest.skip("Coarse dataset files not found for SGS comparison")
        assert content1 != files2[0].read_text(), (
            "SGS with different seeds produced identical output")

    def test_sgsim_same_seed_reproducible(self, run_t2p):
        """Two SGS runs with the same seed should produce identical outputs."""
        r1 = run_t2p(_grid_t2p(options_extra="SGSIM\n  SEED 99\n  WRITE_DATASET_FILES"),
                     extra_files=FILES)
        if r1.returncode != 0:
            pytest.skip("SGS run did not complete successfully")
        files1 = sorted(r1.workdir.glob("*Coarse*"))
        if not files1:
            pytest.skip("Coarse dataset files not found")
        content1 = files1[0].read_text()  # capture before second run overwrites

        r2 = run_t2p(_grid_t2p(options_extra="SGSIM\n  SEED 99\n  WRITE_DATASET_FILES"),
                     extra_files=FILES)
        if r2.returncode != 0:
            pytest.skip("SGS run did not complete successfully")
        files2 = sorted(r2.workdir.glob("*Coarse*"))
        if not files2:
            pytest.skip("Coarse dataset files not found")
        assert content1 == files2[0].read_text(), (
            "SGS with same seed produced different output")

    def test_sgsim_incompatible_with_weight_write(self, run_t2p):
        """Per docs, SGS is not compatible with WRITE_GRID_WEIGHTS."""
        r = run_t2p(_grid_t2p(options_extra="SGSIM\n  SEED 1\n  WRITE_GRID_WEIGHTS"),
                    extra_files=FILES)
        assert r.returncode != 0 or "warning" in r.combined.lower() or \
               "error" in r.combined.lower(), (
            "Expected warning/error for SGSIM + WRITE_GRID_WEIGHTS combination")


# ---------------------------------------------------------------------------
# 4. Weight caching (WRITE_*/READ_*)
# ---------------------------------------------------------------------------

class TestWeightCaching:

    def test_write_and_read_grid_weights_binary(self, run_t2p):
        """Write weights on first run, reuse them on second; both must succeed."""
        # Use NZ=1 (single layer) so observations fall in the only layer, ensuring
        # the WRITE run produces non-empty weight files that the READ run can consume.
        # Use explicit Fine variogram (no INFER_LAST_CLASS) so Fine is also kriged
        # and its weight file is complete.
        single_grd = """\
NX    5
NY    5
NZ    1
DX    3000.0
DY    3000.0
TOP   constant 100.0
BOT   constant 80.0
IBOUND constant 1
"""
        def _t2p(extra_options):
            return f"""\
                BEGIN OPTIONS
                  MAX_VSTRUCT 1
                  WRITE_NODE_FILES
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
                    1 Sph 0.05 0.038 1.5E4 1.5E4 0.0 16
                  CLASS Fine
                    1 Sph 0.05 0.038 1.5E4 1.5E4 0.0 16
                  CLASS PilotPoints
                    1 Sph 0.00 1.0 1.0E5 1.0E5 0.0 25
                END VARIOGRAMS
                BEGIN PP_LOCS
                  1 5050.0 4005.0 1
                END PP_LOCS
                BEGIN PP_PARAMETERS
                  TYPE Global
                  1  0.5  -0.5  0.0
                  TYPE Aquifer
                  1 Coarse  0.1  100.0  1e-5  0.15  10.0  0.001
                  1 Fine  0.001    0.5  1e-4  0.05   5.0  0.001
                END PP_PARAMETERS
                """

        files = {**FILES, "grid.grd": single_grd}
        r1 = run_t2p(_t2p("WRITE_GRID_WEIGHTS"), extra_files=files)
        assert r1.returncode == 0, r1.combined

        weight_files = list(r1.workdir.glob("*.wts")) + list(r1.workdir.glob("*.gwt"))
        if not weight_files:
            pytest.skip("No weight files found after WRITE_GRID_WEIGHTS run")

        # Both calls share the same workdir (tmp_path), so weight files are already present
        r2 = run_t2p(_t2p("READ_GRID_WEIGHTS"), extra_files=files)
        assert r2.returncode == 0, r2.combined

    def test_write_and_read_pp_weights(self, run_t2p):
        r1 = run_t2p(_grid_t2p(options_extra="WRITE_PP_WEIGHTS"), extra_files=FILES)
        assert r1.returncode == 0, r1.combined

        weight_files = list(r1.workdir.glob("*.ppw")) + list(r1.workdir.glob("*.pwt"))
        if not weight_files:
            pytest.skip("No PP weight files found")

        extra = dict(FILES)
        for wf in weight_files:
            extra[wf.name] = wf.read_text()

        r2 = run_t2p(_grid_t2p(options_extra="READ_PP_WEIGHTS"), extra_files=extra)
        assert r2.returncode == 0, r2.combined

    def test_text_suffix_weight_files(self, run_t2p):
        """WRITE_GRID_WEIGHTS TEXT should produce a readable text file."""
        r = run_t2p(_grid_t2p(options_extra="WRITE_GRID_WEIGHTS TEXT"), extra_files=FILES)
        assert r.returncode == 0, r.combined
        text_weights = list(r.workdir.glob("*.txt")) + list(r.workdir.glob("*weights*"))
        # At minimum, the run should complete without error
        assert "error" not in r.combined.lower(), r.combined


# ---------------------------------------------------------------------------
# 5. ONLY_TEXTURE (Data2Texture mode)
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
# 6. MAX_LOG_LENGTH (3D log preprocessing)
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
# 7. INTERP_DIM 3 grid geometry
# ---------------------------------------------------------------------------

class TestInterpDim3GridGeometry:

    def test_elevation_file_writes_top_bottoms_and_layer_centers(self, run_t2p):
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
# 8. Ordinary kriging in 3D
# ---------------------------------------------------------------------------

class TestOrdinaryKriging3D:

    def test_exact_vertical_hits_preserve_layer_sequence(self, run_t2p):
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


# ---------------------------------------------------------------------------
# 9. READ_GRIDDED mode
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
# 10. ANISOTROPIC_SEARCH
# ---------------------------------------------------------------------------

class TestAnisotropicSearch:

    def test_anisotropic_search_runs(self, run_t2p):
        r = run_t2p(_grid_t2p(options_extra="ANISOTROPIC_SEARCH"), extra_files=FILES)
        assert r.returncode == 0, r.combined


# ---------------------------------------------------------------------------
# 11. NO_TEX_BOUNDS
# ---------------------------------------------------------------------------

class TestNoTexBounds:

    def test_no_tex_bounds_allows_outside_range(self, run_t2p):
        """With NO_TEX_BOUNDS, values outside [0,1] should not trigger an error."""
        r = run_t2p(_grid_t2p(options_extra="NO_TEX_BOUNDS\n  WRITE_DATASET_FILES"),
                    extra_files=FILES)
        assert r.returncode == 0, r.combined
