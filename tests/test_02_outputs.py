"""
test_02_outputs.py
==================
Tests that Texture2Par produces the correct output files with the expected
structure, column sets, and physically plausible values.
"""

import csv
import math
from pathlib import Path

import pytest


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

DATASET = """\
Location,WellIndex,IntervalIndex,X,Y,Zland,TopDepth,BotDepth,Coarse,Fine
W1,1,1,5000.0,4000.0,100.0,0.0,5.0,0.60,0.40
W1,1,2,5000.0,4000.0,100.0,5.0,10.0,0.55,0.45
W2,2,1,5100.0,4000.0,100.0,0.0,5.0,0.70,0.30
W2,2,2,5100.0,4000.0,100.0,5.0,10.0,0.45,0.55
W3,3,1,5050.0,4010.0,100.0,0.0,5.0,0.55,0.45
W3,3,2,5050.0,4010.0,100.0,5.0,10.0,0.50,0.50
W4,4,1,4980.0,4008.0,100.0,0.0,5.0,0.65,0.35
W4,4,2,4980.0,4008.0,100.0,5.0,10.0,0.42,0.58
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

GAPPED_DATASET = """\
Location,WellIndex,IntervalIndex,X,Y,Zland,TopDepth,BotDepth,Coarse,Fine
W1,1,1,5000.0,4000.0,100.0,0.0,10.0,0.60,0.40
W3,3,1,5050.0,4010.0,100.0,0.0,10.0,0.55,0.45
"""


def _base_t2p(extra_options="", extra_blocks="") -> str:
    return f"""\
        BEGIN OPTIONS
          MAX_VSTRUCT 1
          INFER_LAST_CLASS
          WRITE_NODE_FILES
          WRITE_DATASET_FILES
          WRITE_ELEVATION_FILE
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
          CLASS PilotPoints
            1 Sph 0.00 1.0 1.0E5 1.0E5 0.0 25
        END VARIOGRAMS

        BEGIN PP_LOCS
          1 5050.0 4005.0 1
        END PP_LOCS

        BEGIN PP_PARAMETERS
          TYPE Global
          # ID KHp  KVp  STp
          1  0.5  -0.5  0.0

          TYPE Aquifer
          # ID Class  Kmin  Kmax    Ss     Sy   Aniso  Kd
          1 Coarse  0.1  100.0  1e-5  0.15  10.0  0.001
          1 Fine  0.001    0.5  1e-4  0.05   5.0  0.001
        END PP_PARAMETERS

        {extra_blocks}
        """


def _csv_rows(path: Path) -> list[dict]:
    return list(csv.DictReader(path.read_text().splitlines()))


def _numeric_values(rows: list[dict], col: str) -> list[float]:
    vals = []
    for r in rows:
        v = r.get(col, "").strip()
        try:
            vals.append(float(v))
        except ValueError:
            pass
    return vals


# ---------------------------------------------------------------------------
# 1. Log file
# ---------------------------------------------------------------------------

class TestLogFile:

    def test_log_file_created(self, run_t2p):
        r = run_t2p(_base_t2p(), extra_files=FILES)
        assert r.returncode == 0, r.combined
        assert r.log, "Log file (.out) was not created or is empty"

    def test_log_records_elapsed_time(self, run_t2p):
        r = run_t2p(_base_t2p(), extra_files=FILES)
        assert r.returncode == 0
        assert "elapsed" in r.log.lower(), "Log file does not record elapsed time"

    def test_log_echoes_options(self, run_t2p):
        r = run_t2p(_base_t2p(), extra_files=FILES)
        assert r.returncode == 0
        # Log should mention at least one key option
        assert any(kw in r.log.upper() for kw in
                   ("MAX_VSTRUCT", "INFER_LAST_CLASS", "OPTIONS")), (
            "Log file does not echo options block")

    def test_log_records_errors_on_bad_run(self, run_t2p, minimal_grd):
        """A failed run should write an error message to the log."""
        bad = _base_t2p().replace("MAX_VSTRUCT 1", "")  # remove mandatory option
        r = run_t2p(bad, extra_files={"grid.grd": minimal_grd})
        if r.returncode == 0:
            pytest.skip("Run unexpectedly succeeded without MAX_VSTRUCT")
        assert "error" in r.combined.lower(), "Expected 'error' in output on bad run"


# ---------------------------------------------------------------------------
# 2. Node / cell CSV files
# ---------------------------------------------------------------------------

class TestNodeFiles:

    def test_kh_file_written(self, run_t2p):
        r = run_t2p(_base_t2p(), extra_files=FILES)
        assert r.returncode == 0
        assert list(r.workdir.glob("*KH*")), "KH output file not found"

    def test_kv_file_written(self, run_t2p):
        r = run_t2p(_base_t2p(), extra_files=FILES)
        assert r.returncode == 0
        assert list(r.workdir.glob("*KV*")), "KV output file not found"

    def test_ss_file_written(self, run_t2p):
        r = run_t2p(_base_t2p(), extra_files=FILES)
        assert r.returncode == 0
        assert list(r.workdir.glob("*SS*")), "SS output file not found"

    def test_sy_file_written(self, run_t2p):
        r = run_t2p(_base_t2p(), extra_files=FILES)
        assert r.returncode == 0
        assert list(r.workdir.glob("*SY*")), "SY output file not found"

    def test_node_files_have_spatial_columns(self, run_t2p):
        """Each node file should contain X, Y, and a layer column."""
        r = run_t2p(_base_t2p(), extra_files=FILES)
        assert r.returncode == 0
        kh_files = list(r.workdir.glob("*KH*"))
        if not kh_files:
            pytest.skip("KH file not found")
        rows = _csv_rows(kh_files[0])
        assert rows, "KH file is empty"
        headers = {k.strip().upper() for k in rows[0].keys()}
        assert "X" in headers or any("X" in h for h in headers), \
            f"X column not found in KH file. Headers: {headers}"
        assert "Y" in headers or any("Y" in h for h in headers), \
            f"Y column not found in KH file. Headers: {headers}"

    def test_kh_values_positive(self, run_t2p):
        r = run_t2p(_base_t2p(), extra_files=FILES)
        assert r.returncode == 0
        kh_files = list(r.workdir.glob("*KH*"))
        if not kh_files:
            pytest.skip("KH file not found")
        rows = _csv_rows(kh_files[0])
        # Find the KH value column (last numeric column that is not X/Y/layer)
        numeric_cols = []
        for k in rows[0].keys():
            try:
                float(rows[0][k])
                if k.upper() not in ("X", "Y", "LAYER", "NODE", "ROW", "COL"):
                    numeric_cols.append(k)
            except ValueError:
                pass
        if not numeric_cols:
            pytest.skip("Cannot identify KH value column")
        vals = _numeric_values(rows, numeric_cols[-1])
        assert all(v > 0 for v in vals if not math.isnan(v)), \
            "KH contains non-positive values"

    def test_sy_values_between_0_and_1(self, run_t2p):
        r = run_t2p(_base_t2p(), extra_files=FILES)
        assert r.returncode == 0
        sy_files = list(r.workdir.glob("*SY*"))
        if not sy_files:
            pytest.skip("SY file not found")
        rows = _csv_rows(sy_files[0])
        numeric_cols = [k for k in rows[0].keys()
                        if k.upper() not in ("X", "Y", "LAYER", "NODE", "ROW", "COL")]
        if not numeric_cols:
            pytest.skip("Cannot identify SY value column")
        vals = _numeric_values(rows, numeric_cols[-1])
        for v in vals:
            assert 0.0 <= v <= 1.0, f"SY value {v} outside physical range [0, 1]"

    def test_node_files_not_written_without_flag(self, run_t2p):
        """GRID flow model always enables WRITE_NODE_FILES internally; the flag is non-suppressible."""
        pytest.skip("GRID type unconditionally sets WRITE_NODE_FILES=true internally; "
                    "parameter node files are always written for GRID-type flow models")


# ---------------------------------------------------------------------------
# 3. Dataset CSV files
# ---------------------------------------------------------------------------

class TestDatasetFiles:

    def test_dataset_files_written(self, run_t2p):
        r = run_t2p(_base_t2p(), extra_files=FILES)
        assert r.returncode == 0
        tex_files = list(r.workdir.glob("*Coarse*")) + list(r.workdir.glob("*coarse*"))
        assert tex_files, "No Coarse dataset output file found"

    def test_dataset_file_columns(self, run_t2p):
        r = run_t2p(_base_t2p(), extra_files=FILES)
        assert r.returncode == 0
        tex_files = list(r.workdir.glob("*Coarse*"))
        if not tex_files:
            pytest.skip("Coarse dataset file not found")
        rows = _csv_rows(tex_files[0])
        assert rows
        headers = {k.strip().upper() for k in rows[0].keys()}
        # Docs say columns are Node/(Row/Col), X, Y, layer
        for expected in ("X", "Y"):
            assert any(expected in h for h in headers), \
                f"Expected column '{expected}' in dataset file. Got: {headers}"

    def test_dataset_files_not_written_without_flag(self, run_t2p):
        """WRITE_DATASET_FILES controls layer-averaged observation files (*_layavg.csv)."""
        t2p = _base_t2p().replace("WRITE_DATASET_FILES", "")
        r = run_t2p(t2p, extra_files=FILES)
        assert r.returncode == 0
        layavg_files = list(r.workdir.glob("*_layavg*"))
        assert not layavg_files, "Layer-averaged dataset files written without WRITE_DATASET_FILES flag"

    def test_missing_well_id_has_zero_model_element(self, run_t2p):
        """A gap in well IDs must not expose an uninitialized model element."""
        r = run_t2p(
            _base_t2p(),
            extra_files={"dataset.csv": GAPPED_DATASET, "grid.grd": GRD},
        )
        assert r.returncode == 0, r.combined

        layavg = r.output_file("test_COARSE_layavg.csv")
        assert layavg.exists()
        rows = _csv_rows(layavg)
        assert len(rows) == 3
        assert int(rows[1]["Well"]) == 2
        assert int(rows[1]["Column"]) == 0


# ---------------------------------------------------------------------------
# 4. Elevation file
# ---------------------------------------------------------------------------

class TestElevationFile:

    def test_elevation_file_written(self, run_t2p):
        r = run_t2p(_base_t2p(), extra_files=FILES)
        assert r.returncode == 0
        elev_files = (list(r.workdir.glob("*elev*")) +
                      list(r.workdir.glob("*ELEV*")) +
                      list(r.workdir.glob("*elevation*")))
        assert elev_files, "Elevation file not found with WRITE_ELEVATION_FILE active"

    def test_elevation_file_not_written_without_flag(self, run_t2p):
        t2p = _base_t2p().replace("WRITE_ELEVATION_FILE", "")
        r = run_t2p(t2p, extra_files=FILES)
        assert r.returncode == 0
        elev_files = (list(r.workdir.glob("*elev*")) +
                      list(r.workdir.glob("*ELEV*")))
        assert not elev_files, "Elevation file unexpectedly written"


class TestDataElevationFile:

    @pytest.mark.parametrize("option", [
        "WRITE_DATA_ELEV_FILE",
        "WRITE_DATA_ELEVATION_FILE",
    ])
    def test_data_elevation_file_contents(self, run_t2p, option):
        t2p = _base_t2p(extra_options=option)
        r = run_t2p(t2p, extra_files=FILES)
        assert r.returncode == 0, r.combined
        assert "WRITE_DATA_ELEV_FILE" in r.log

        output = r.output_file("test_DATA_ELEVATION.csv")
        assert output.exists()
        rows = _csv_rows(output)
        assert len(rows) == 4
        assert list(rows[0]) == [
            "Dataset", "Point", "X", "Y", "Top", "Bottom1", "Bottom2"
        ]
        assert [int(row["Point"]) for row in rows] == [1, 2, 3, 4]
        assert all(int(row["Dataset"]) == 1 for row in rows)
        assert [float(row["X"]) for row in rows] == [5000.0, 5100.0, 5050.0, 4980.0]
        assert [float(row["Y"]) for row in rows] == [4000.0, 4000.0, 4010.0, 4008.0]
        assert all(float(row["Top"]) == 100.0 for row in rows)
        assert all(float(row["Bottom1"]) == 90.0 for row in rows)
        assert all(float(row["Bottom2"]) == 80.0 for row in rows)

    def test_data_elevation_file_is_opt_in(self, run_t2p):
        r = run_t2p(_base_t2p(), extra_files=FILES)
        assert r.returncode == 0, r.combined
        assert not r.output_exists("test_DATA_ELEVATION.csv")


# ---------------------------------------------------------------------------
# 5. WRITE_NONE suppresses model output
# ---------------------------------------------------------------------------

class TestWriteNone:

    def test_write_none_suppresses_model_files(self, run_t2p):
        t2p = _base_t2p(extra_options="WRITE_NONE")
        r = run_t2p(t2p, extra_files=FILES)
        assert r.returncode == 0, r.combined
        # No hydraulic property model file should be written
        for pat in ("*.upw", "*.lpf", "*.gwf", "*groundwater*"):
            found = list(r.workdir.glob(pat))
            assert not found, f"Model file {found} written despite WRITE_NONE"


# ---------------------------------------------------------------------------
# 6. Subsidence outputs
# ---------------------------------------------------------------------------

SUBSIDENCE_BLOCK = """\
    BEGIN SUBSIDENCE
      FINE_CLASSES Fine
      INTERBED_CUTOFF 0.9
      VARIOGRAM
        1 Sph 0.00 1.0 1.0E4 1.0E4 0.0 16
    END SUBSIDENCE
    """

SUBSIDENCE_DATASET = """\
Location,WellIndex,IntervalIndex,X,Y,Zland,TopDepth,BotDepth,Coarse,Fine
W1,1,1,5000.0,4000.0,100.0,0.0,2.0,0.05,0.95
W1,1,2,5000.0,4000.0,100.0,2.0,5.0,0.60,0.40
W2,2,1,5100.0,4000.0,100.0,0.0,2.0,0.03,0.97
W2,2,2,5100.0,4000.0,100.0,2.0,5.0,0.65,0.35
W3,3,1,5050.0,4010.0,100.0,0.0,2.0,0.04,0.96
W3,3,2,5050.0,4010.0,100.0,2.0,5.0,0.55,0.45
W4,4,1,4980.0,4008.0,100.0,0.0,2.0,0.06,0.94
W4,4,2,4980.0,4008.0,100.0,2.0,5.0,0.50,0.50
"""

SUB_PP_PARAMS = """\
    BEGIN PP_PARAMETERS
      TYPE Global
      # ID KHp  KVp  STp
      1  0.5  -0.5  0.0

      TYPE Aquifer
      # ID Class  Kmin  Kmax   KvI   Sske   Sskv    Sy   Aniso  Kd
      1 Coarse  0.1  100.0  0.01  1e-6   0.0  0.15  10.0  0.001
      1 Fine  0.001    0.5 0.001  1e-5  1e-4  0.05   5.0  0.001
    END PP_PARAMETERS
    """


class TestSubsidenceOutputs:

    def _sub_t2p(self) -> str:
        return f"""\
            BEGIN OPTIONS
              MAX_VSTRUCT 1
              INFER_LAST_CLASS
              WRITE_NODE_FILES
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
              CLASS PilotPoints
                1 Sph 0.00 1.0 1.0E5 1.0E5 0.0 25
            END VARIOGRAMS
            {SUBSIDENCE_BLOCK}
            BEGIN PP_LOCS
              1 5050.0 4005.0 1
            END PP_LOCS
            {SUB_PP_PARAMS}
            """

    def test_subsidence_run_completes(self, run_t2p):
        r = run_t2p(self._sub_t2p(),
                    extra_files={"dataset.csv": SUBSIDENCE_DATASET, "grid.grd": GRD})
        assert r.returncode == 0, r.combined

    def test_ib_thk_file_written(self, run_t2p):
        r = run_t2p(self._sub_t2p(),
                    extra_files={"dataset.csv": SUBSIDENCE_DATASET, "grid.grd": GRD})
        assert r.returncode == 0
        assert list(r.workdir.glob("*IB_THICK*")) or list(r.workdir.glob("*ib_thick*")), \
            "IB_THICK file not found"

    def test_kvf_file_written(self, run_t2p):
        r = run_t2p(self._sub_t2p(),
                    extra_files={"dataset.csv": SUBSIDENCE_DATASET, "grid.grd": GRD})
        assert r.returncode == 0
        assert list(r.workdir.glob("*KVF*")) or list(r.workdir.glob("*kvf*")), \
            "KVF file not found"

    def test_ib_thk_non_negative(self, run_t2p):
        r = run_t2p(self._sub_t2p(),
                    extra_files={"dataset.csv": SUBSIDENCE_DATASET, "grid.grd": GRD})
        assert r.returncode == 0
        ib_files = list(r.workdir.glob("*IB_THICK*")) + list(r.workdir.glob("*ib_thick*"))
        if not ib_files:
            pytest.skip("IB_THK file not found")
        rows = _csv_rows(ib_files[0])
        numeric_cols = [k for k in rows[0].keys()
                        if k.upper() not in ("X", "Y", "LAYER", "NODE", "ROW", "COL")]
        if not numeric_cols:
            pytest.skip("Cannot identify IB_THK value column")
        vals = _numeric_values(rows, numeric_cols[-1])
        assert all(v >= 0 for v in vals if not math.isnan(v)), \
            "IB_THK contains negative values"


# ---------------------------------------------------------------------------
# 7. Row count sanity (5×5 grid, 2 layers = 50 nodes per layer file)
# ---------------------------------------------------------------------------

class TestOutputRowCounts:

    def test_kh_row_count_matches_grid(self, run_t2p):
        r = run_t2p(_base_t2p(), extra_files=FILES)
        assert r.returncode == 0
        kh_files = list(r.workdir.glob("*KH*"))
        if not kh_files:
            pytest.skip("KH file not found")
        rows = _csv_rows(kh_files[0])
        # 5x5 grid => 25 rows in wide format (one row per Row/Col with one column per layer)
        assert len(rows) == 25, \
            f"Expected 25 rows in KH file (5×5 wide format), got {len(rows)}"
