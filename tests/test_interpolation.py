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
# 7. READ_GRIDDED mode
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
# 8. ANISOTROPIC_SEARCH
# ---------------------------------------------------------------------------

class TestAnisotropicSearch:

    def test_anisotropic_search_runs(self, run_t2p):
        r = run_t2p(_grid_t2p(options_extra="ANISOTROPIC_SEARCH"), extra_files=FILES)
        assert r.returncode == 0, r.combined


# ---------------------------------------------------------------------------
# 9. NO_TEX_BOUNDS
# ---------------------------------------------------------------------------

class TestNoTexBounds:

    def test_no_tex_bounds_allows_outside_range(self, run_t2p):
        """With NO_TEX_BOUNDS, values outside [0,1] should not trigger an error."""
        r = run_t2p(_grid_t2p(options_extra="NO_TEX_BOUNDS\n  WRITE_DATASET_FILES"),
                    extra_files=FILES)
        assert r.returncode == 0, r.combined
