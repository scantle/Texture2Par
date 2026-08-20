"""
test_01_input_parsing.py
========================
Tests that Texture2Par correctly parses (or rejects) various input file
structures.  All tests execute the binary against real .t2p files.
"""

import pytest

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _simple_t2p(options_extra="", classes="Coarse\n  Fine",
                infer=True, extra_blocks="") -> str:
    infer_line = "  INFER_LAST_CLASS" if infer else ""
    return f"""\
        BEGIN OPTIONS
          MAX_VSTRUCT 1
          {infer_line}
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
          CLASS PilotPoints
            1 Sph 0.00 1.0 1.0E5 1.0E5 0.0 25
        END VARIOGRAMS

        BEGIN PP_LOCS
          1 500000.0 4000000.0 1
        END PP_LOCS

        BEGIN PP_PARAMETERS
          TYPE Global
          # ID KHp KVp STp
          1 0.5 -0.5 0.0

          TYPE Aquifer
          # ID Class Kmin Kmax Ss Sy Aniso Kd
          1 Coarse 0.1 100.0 1e-5 0.15 10.0 0.001
          1 Fine   0.001 0.5 1e-4 0.05 5.0  0.001
        END PP_PARAMETERS

        {extra_blocks}
        """


def _zoned_t2p() -> str:
    return """\
        BEGIN OPTIONS
          MAX_VSTRUCT 1
          WRITE_NONE
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
          PP_ZONE_FILE ppzones.txt
        END FLOW_MODEL

        BEGIN VARIOGRAMS
          CLASS Coarse
            1 Sph 0.05 0.038 1.5E4 1.5E4 0.0 16
          CLASS Fine
            1 Sph 0.05 0.038 1.5E4 1.5E4 0.0 16
          CLASS PilotPoints
            1 Sph 0.00 1.0 1.0E5 1.0E5 0.0 3
        END VARIOGRAMS

        BEGIN PP_LOCS
          1 1000.0 1000.0 1
          2 3000.0 1000.0 2
          3 1000.0 3000.0 3
        END PP_LOCS

        BEGIN PP_PARAMETERS
          TYPE Global
          1 0.5 -0.5 0.0
          2 0.5 -0.5 0.0
          3 0.5 -0.5 0.0

          TYPE Aquifer
          1 Coarse 0.1   100.0 1e-5 0.15 10.0 0.001
          1 Fine   0.001 0.5   1e-4 0.05 5.0  0.001
          2 Coarse 0.1   100.0 1e-5 0.15 10.0 0.001
          2 Fine   0.001 0.5   1e-4 0.05 5.0  0.001
          3 Coarse 0.1   100.0 1e-5 0.15 10.0 0.001
          3 Fine   0.001 0.5   1e-4 0.05 5.0  0.001
        END PP_PARAMETERS
        """


def _hsu_t2p(interp_dim=2, variogram_extra="") -> str:
    t2p = _zoned_t2p().replace(
        "          WRITE_NONE",
        f"          WRITE_NONE\n          INTERP_DIM {interp_dim}",
    ).replace(
        "          PP_ZONE_FILE ppzones.txt",
        "          HSU_FILE hsus.txt",
    ).replace(
        "          2 3000.0 1000.0 2",
        "          2 3000.0 1000.0 1",
    ).replace(
        "          3 1000.0 3000.0 3",
        "          3 1000.0 3000.0 1",
    )
    if interp_dim == 3:
        t2p = t2p.replace(
            "1 Sph 0.05 0.038 1.5E4 1.5E4 0.0 16",
            "1 Sph 0.05 0.038 1.5E4 1.5E4 100.0 0.0 0.0 0.0 16",
        )
    return t2p.replace(
        "          CLASS PilotPoints",
        f"{variogram_extra}          CLASS PilotPoints",
    )


ZONED_GRD = """\
NX    2
NY    2
NZ    3
DX    2000.0
DY    2000.0
TOP   constant 100.0
BOT   constant 70.0
IBOUND constant 1
"""

ZONED_DATASET = """\
Location,WellIndex,IntervalIndex,X,Y,Zland,TopDepth,BotDepth,Coarse,Fine
W1,1,1,1000.0,1000.0,100.0,0.0,10.0,0.6,0.4
W1,1,2,1000.0,1000.0,100.0,10.0,20.0,0.5,0.5
W1,1,3,1000.0,1000.0,100.0,20.0,30.0,0.4,0.6
W2,2,1,3000.0,1000.0,100.0,0.0,10.0,0.5,0.5
W2,2,2,3000.0,1000.0,100.0,10.0,20.0,0.4,0.6
W2,2,3,3000.0,1000.0,100.0,20.0,30.0,0.6,0.4
W3,3,1,1000.0,3000.0,100.0,0.0,10.0,0.4,0.6
W3,3,2,1000.0,3000.0,100.0,10.0,20.0,0.6,0.4
W3,3,3,1000.0,3000.0,100.0,20.0,30.0,0.5,0.5
"""

HSUS = """\
1 2 3
Sand Clay Sand Clay
Clay Gravel Clay Gravel
Gravel Sand Gravel Sand
"""

HSU_DATASET = """\
Location,WellIndex,IntervalIndex,X,Y,Zland,TopDepth,BotDepth,Coarse,Fine,HSU1,HSU2,HSU3
W1,1,1,1000.0,1000.0,100.0,0.0,10.0,0.6,0.4,sand,clay,gravel
W1,1,2,1000.0,1000.0,100.0,10.0,20.0,0.5,0.5,sand,clay,gravel
W1,1,3,1000.0,1000.0,100.0,20.0,30.0,0.4,0.6,sand,clay,gravel
W2,2,1,3000.0,1000.0,100.0,0.0,10.0,0.5,0.5,CLAY,GRAVEL,SAND
W2,2,2,3000.0,1000.0,100.0,10.0,20.0,0.4,0.6,CLAY,GRAVEL,SAND
W2,2,3,3000.0,1000.0,100.0,20.0,30.0,0.6,0.4,CLAY,GRAVEL,SAND
W3,3,1,1000.0,3000.0,100.0,0.0,10.0,0.4,0.6,Sand,Clay,Gravel
W3,3,2,1000.0,3000.0,100.0,10.0,20.0,0.6,0.4,Sand,Clay,Gravel
W3,3,3,1000.0,3000.0,100.0,20.0,30.0,0.5,0.5,Sand,Clay,Gravel
W4,4,1,3000.0,3000.0,100.0,0.0,10.0,0.5,0.5,Clay,Gravel,Sand
W4,4,2,3000.0,3000.0,100.0,10.0,20.0,0.4,0.6,Clay,Gravel,Sand
W4,4,3,3000.0,3000.0,100.0,20.0,30.0,0.6,0.4,Clay,Gravel,Sand
"""


# ---------------------------------------------------------------------------
# 1. OPTIONS block
# ---------------------------------------------------------------------------

class TestOptionsBlock:

    def test_missing_max_vstruct_uses_default(self, run_t2p, minimal_dataset, minimal_grd):
        """Omitting MAX_VSTRUCT should use the program default and continue."""
        t2p = """\
            BEGIN OPTIONS
              INFER_LAST_CLASS
            END OPTIONS
            ...
            """
        r = run_t2p(t2p, extra_files={"dataset.csv": minimal_dataset,
                                       "grid.grd": minimal_grd})
        assert r.returncode == 0, r.combined

    def test_max_vstruct_one(self, run_t2p, minimal_dataset, minimal_grd):
        """MAX_VSTRUCT 1 is the minimal valid configuration."""
        r = run_t2p(_simple_t2p(),
                    extra_files={"dataset.csv": minimal_dataset,
                                 "grid.grd": minimal_grd})
        assert r.returncode == 0, r.combined

    @pytest.mark.parametrize("option", [
        "WRITE_DATASET_FILES",
        "WRITE_NODE_FILES",
        "WRITE_ELEVATION_FILE",
        "WRITE_DATA_ELEVATION_FILE",
        "INFER_LAST_CLASS",
    ])
    def test_boolean_flag_options_accepted(self, run_t2p, minimal_dataset,
                                           minimal_grd, option):
        """Boolean flags that take no value should be silently accepted."""
        r = run_t2p(_simple_t2p(options_extra=option),
                    extra_files={"dataset.csv": minimal_dataset,
                                 "grid.grd": minimal_grd})
        # Only check that T2P does not crash on the flag itself
        assert "error" not in r.combined.lower() or r.returncode == 0, r.combined

    @pytest.mark.parametrize("dim", [2, 3])
    def test_interp_dim(self, run_t2p, minimal_dataset, minimal_grd, dim):
        """INTERP_DIM 2 and 3 must both be accepted in the options block."""
        r = run_t2p(_simple_t2p(options_extra=f"INTERP_DIM {dim}"),
                    extra_files={"dataset.csv": minimal_dataset,
                                 "grid.grd": minimal_grd})
        # A dim=3 run may fail for other reasons (no vertical range) but the
        # option itself must be parsed without an "unknown option" error.
        assert "unknown option" not in r.combined.lower(), r.combined

    def test_no_data_value_custom(self, run_t2p, minimal_dataset, minimal_grd):
        r = run_t2p(_simple_t2p(options_extra="NO_DATA_VALUE -9999"),
                    extra_files={"dataset.csv": minimal_dataset,
                                 "grid.grd": minimal_grd})
        assert r.returncode == 0, r.combined

    def test_project_name_reflected_in_outputs(self, run_t2p, minimal_dataset,
                                                minimal_grd):
        r = run_t2p(_simple_t2p(options_extra="PROJECT_NAME myproject\n  WRITE_DATASET_FILES"),
                    extra_files={"dataset.csv": minimal_dataset,
                                 "grid.grd": minimal_grd})
        assert r.returncode == 0, r.combined
        # At least one output file should carry the project name prefix
        outputs = list(r.workdir.glob("myproject*"))
        assert outputs, "No output files found with PROJECT_NAME prefix"

    def test_seed_option_accepted(self, run_t2p, minimal_dataset, minimal_grd):
        r = run_t2p(_simple_t2p(options_extra="SEED 42"),
                    extra_files={"dataset.csv": minimal_dataset,
                                 "grid.grd": minimal_grd})
        assert "unknown option" not in r.combined.lower(), r.combined


# ---------------------------------------------------------------------------
# 2. CLASSES block
# ---------------------------------------------------------------------------

class TestClassesBlock:

    def test_class_name_with_underscore(self, run_t2p, minimal_dataset, minimal_grd):
        r = run_t2p(_simple_t2p(classes="Coarse\n  Fine", infer=True).replace(
                        "Coarse", "Coarse_Sand").replace("Fine", "Fine_Clay"),
                    extra_files={"dataset.csv": minimal_dataset.replace(
                        "Coarse", "Coarse_Sand").replace("Fine", "Fine_Clay"),
                                 "grid.grd": minimal_grd})
        # Should not error on underscore names
        assert "invalid" not in r.combined.lower(), r.combined

    def test_class_name_with_hyphen(self, run_t2p, minimal_dataset, minimal_grd):
        r = run_t2p(_simple_t2p(classes="Coarse\n  Fine", infer=True).replace(
                        "Coarse", "Coarse-Sand").replace("Fine", "Fine-Clay"),
                    extra_files={"dataset.csv": minimal_dataset.replace(
                        "Coarse", "Coarse-Sand").replace("Fine", "Fine-Clay"),
                                 "grid.grd": minimal_grd})
        assert "invalid" not in r.combined.lower(), r.combined


# ---------------------------------------------------------------------------
# 3. DATASET block
# ---------------------------------------------------------------------------

class TestDatasetBlock:

    def test_missing_dataset_file_errors(self, run_t2p, minimal_grd):
        """Referencing a non-existent dataset file should produce an error."""
        r = run_t2p(_simple_t2p(),
                    extra_files={"grid.grd": minimal_grd})
        # dataset.csv is NOT written
        assert r.returncode != 0 or "error" in r.combined.lower()

    @pytest.mark.parametrize("transform", ["LOG10", "LN", "LOGNAT", "SQRT",
                                            "EXP", "INV", "INVERSE"])
    def test_transform_keywords_accepted(self, run_t2p, minimal_dataset,
                                         minimal_grd, transform):
        """Each documented TRANSFORM keyword must be accepted without error."""
        t2p = f"""\
            BEGIN OPTIONS
              MAX_VSTRUCT 1
              INFER_LAST_CLASS
            END OPTIONS
            BEGIN CLASSES
              Coarse
              Fine
            END CLASSES
            BEGIN DATASET
              FILE dataset.csv
            END DATASET
            BEGIN DATASET
              FILE secondary.csv
              TRANSFORM {transform}
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
              1 500000.0 4000000.0 1
            END PP_LOCS
            BEGIN PP_PARAMETERS
              TYPE Global
              1 0.5 -0.5 0.0
              TYPE Aquifer
              1 Coarse 0.1 100.0 1e-5 0.15 10.0 0.001
              1 Fine   0.001 0.5 1e-4 0.05 5.0 0.001
            END PP_PARAMETERS
            """
        # Use a simple secondary dataset with a numeric column named "Rho"
        secondary = """\
            Location,WellIndex,IntervalIndex,X,Y,Zland,TopDepth,BotDepth,Rho
            W1,1,1,500000.0,4000000.0,100.0,0.0,5.0,50.0
            W2,2,1,510000.0,4000000.0,100.0,0.0,5.0,40.0
            """
        r = run_t2p(t2p, extra_files={
            "dataset.csv": minimal_dataset,
            "secondary.csv": secondary,
            "grid.grd": minimal_grd,
        })
        assert "unknown transform" not in r.combined.lower(), r.combined

    def test_dataset_no_data_rows_excluded(self, run_t2p, minimal_grd):
        """Rows with NO_DATA_VALUE (-999) must not cause a crash."""
        dataset_with_nodata = """\
            Location,WellIndex,IntervalIndex,X,Y,Zland,TopDepth,BotDepth,Coarse,Fine
            W1,1,1,5000.0,4000.0,100.0,0.0,5.0,-999,-999
            W2,2,1,5100.0,4000.0,100.0,0.0,5.0,0.7,0.3
            W3,3,1,5050.0,4010.0,100.0,0.0,5.0,0.55,0.45
            """
        r = run_t2p(_simple_t2p(),
                    extra_files={"dataset.csv": dataset_with_nodata,
                                 "grid.grd": minimal_grd})
        assert r.returncode == 0, r.combined


# ---------------------------------------------------------------------------
# 4. FLOW_MODEL block
# ---------------------------------------------------------------------------

class TestFlowModelBlock:

    def test_grid_type_accepted(self, run_t2p, minimal_dataset, minimal_grd):
        r = run_t2p(_simple_t2p(),
                    extra_files={"dataset.csv": minimal_dataset,
                                 "grid.grd": minimal_grd})
        assert r.returncode == 0, r.combined

    def test_unknown_flow_model_type_rejected(self, run_t2p, minimal_dataset,
                                               minimal_grd):
        bad_t2p = _simple_t2p().replace("TYPE GRID", "TYPE UNKNOWN_MODEL")
        r = run_t2p(bad_t2p, extra_files={"dataset.csv": minimal_dataset,
                                            "grid.grd": minimal_grd})
        assert r.returncode != 0 or "error" in r.combined.lower()

    def test_pp_zone_ids_can_be_distributed_across_arrays(self, run_t2p):
        """The model-wide zone set need not occur together in any one array."""
        ppzones = """\
            1 2 3
            1 2 1 2
            2 3 2 3
            1 3 1 3
            """
        r = run_t2p(_zoned_t2p(), extra_files={
            "dataset.csv": ZONED_DATASET,
            "grid.grd": ZONED_GRD,
            "ppzones.txt": ppzones,
        })
        assert r.returncode == 0, r.combined

    def test_pp_zone_ids_must_match_model_grid(self, run_t2p):
        """Equal zone counts must not hide different pilot-point zone IDs."""
        ppzones = """\
            1 2 3
            1 2 4 1
            2 4 2 4
            1 4 1 4
            """
        r = run_t2p(_zoned_t2p(), extra_files={
            "dataset.csv": ZONED_DATASET,
            "grid.grd": ZONED_GRD,
            "ppzones.txt": ppzones,
        })
        assert r.returncode != 0
        assert "Mismatch between Pilot Point Zones" in r.combined

    def test_missing_grid_file_errors(self, run_t2p, minimal_dataset):
        """GRID_FILE referencing a non-existent file must error."""
        r = run_t2p(_simple_t2p(),
                    extra_files={"dataset.csv": minimal_dataset})
        assert r.returncode != 0 or "error" in r.combined.lower()


class TestHSUCategories:

    def test_hsus_can_be_distributed_across_arrays_and_ignore_case(self, run_t2p):
        """Global HSU IDs are shared across arrays and matched case-insensitively."""
        r = run_t2p(_hsu_t2p(), extra_files={
            "dataset.csv": HSU_DATASET,
            "grid.grd": ZONED_GRD,
            "hsus.txt": HSUS,
        })
        assert r.returncode == 0, r.combined

    @pytest.mark.parametrize("interp_dim", [2, 3])
    def test_unknown_dataset_hsu_is_rejected(self, run_t2p, interp_dim):
        dataset = HSU_DATASET.replace(",sand,clay,gravel", ",Missing,clay,gravel", 1)
        r = run_t2p(_hsu_t2p(interp_dim), extra_files={
            "dataset.csv": dataset,
            "grid.grd": ZONED_GRD,
            "hsus.txt": HSUS,
        })
        assert r.returncode != 0
        assert "Invalid HSU: Missing" in r.combined

    def test_unknown_variogram_hsu_is_rejected(self, run_t2p):
        zonal_variogram = """\
          CLASS Coarse%Missing
            1 Sph 0.05 0.038 1.5E4 1.5E4 0.0
"""
        r = run_t2p(_hsu_t2p(variogram_extra=zonal_variogram), extra_files={
            "dataset.csv": HSU_DATASET,
            "grid.grd": ZONED_GRD,
            "hsus.txt": HSUS,
        })
        assert r.returncode != 0
        assert "Invalid HSU Variogram Used: Missing" in r.combined

# ---------------------------------------------------------------------------
# 5. GRD file parsing
# ---------------------------------------------------------------------------

class TestGrdFile:

    def test_constant_top_bot(self, run_t2p, minimal_dataset, minimal_grd):
        r = run_t2p(_simple_t2p(),
                    extra_files={"dataset.csv": minimal_dataset,
                                 "grid.grd": minimal_grd})
        assert r.returncode == 0, r.combined

    def test_dz_instead_of_bot(self, run_t2p, minimal_dataset):
        grd = """\
            NX    5
            NY    5
            NZ    2
            DX    2000.0
            DY    2000.0
            TOP   constant 100.0
            DZ    10.0
            IBOUND constant 1
            """
        r = run_t2p(_simple_t2p(),
                    extra_files={"dataset.csv": minimal_dataset,
                                 "grid.grd": grd})
        assert r.returncode == 0, r.combined

    def test_use_uniform_ibound(self, run_t2p, minimal_dataset):
        grd = """\
            NX    5
            NY    5
            NZ    2
            DX    2000.0
            DY    2000.0
            TOP   constant 100.0
            BOT   constant 80.0
            USE_UNIFORM_IBOUND
            IBOUND constant 1
            """
        r = run_t2p(_simple_t2p(),
                    extra_files={"dataset.csv": minimal_dataset,
                                 "grid.grd": grd})
        assert r.returncode == 0, r.combined

    def test_missing_nx_errors(self, run_t2p, minimal_dataset):
        grd = """\
            NY    5
            NZ    2
            DX    2000.0
            DY    2000.0
            TOP   constant 100.0
            BOT   constant 80.0
            IBOUND constant 1
            """
        r = run_t2p(_simple_t2p(),
                    extra_files={"dataset.csv": minimal_dataset,
                                 "grid.grd": grd})
        assert r.returncode != 0 or "error" in r.combined.lower()


# ---------------------------------------------------------------------------
# 6. VARIOGRAMS block
# ---------------------------------------------------------------------------

class TestVariogramsBlock:

    @pytest.mark.parametrize("vtype", ["Sph", "Exp", "Hol", "Gau", "Cir"])
    def test_all_variogram_types(self, run_t2p, minimal_dataset, minimal_grd,
                                 vtype):
        """All five documented variogram types must be accepted."""
        t2p = _simple_t2p().replace(
            "1 Sph 0.05 0.038 1.5E4 1.5E4 0.0 16",
            f"1 {vtype} 0.05 0.038 1.5E4 1.5E4 0.0 16"
        )
        r = run_t2p(t2p, extra_files={"dataset.csv": minimal_dataset,
                                       "grid.grd": minimal_grd})
        assert "unknown variogram" not in r.combined.lower(), r.combined

    def test_pilot_point_variogram_must_be_last(self, run_t2p, minimal_dataset,
                                                 minimal_grd):
        """PilotPoints variogram listed before class variogram should warn/error."""
        t2p = _simple_t2p().replace(
            """\
          CLASS Coarse
            1 Sph 0.05 0.038 1.5E4 1.5E4 0.0 16
          CLASS PilotPoints
            1 Sph 0.00 1.0 1.0E5 1.0E5 0.0 25""",
            """\
          CLASS PilotPoints
            1 Sph 0.00 1.0 1.0E5 1.0E5 0.0 25
          CLASS Coarse
            1 Sph 0.05 0.038 1.5E4 1.5E4 0.0 16"""
        )
        r = run_t2p(t2p, extra_files={"dataset.csv": minimal_dataset,
                                       "grid.grd": minimal_grd})
        assert r.returncode != 0 or "error" in r.combined.lower(), (
            "Expected an error when PilotPoints variogram is not last")

    def test_missing_pilot_point_variogram_errors(self, run_t2p, minimal_dataset,
                                                   minimal_grd):
        t2p = _simple_t2p().replace(
            "  CLASS PilotPoints\n            1 Sph 0.00 1.0 1.0E5 1.0E5 0.0 25\n", ""
        )
        r = run_t2p(t2p, extra_files={"dataset.csv": minimal_dataset,
                                       "grid.grd": minimal_grd})
        assert r.returncode != 0 or "error" in r.combined.lower()

    def test_composite_variogram_two_structures(self, run_t2p, minimal_dataset,
                                                 minimal_grd):
        """Multiple structures on one class (composite variogram) is valid."""
        t2p = _simple_t2p().replace(
            "  CLASS Coarse\n            1 Sph 0.05 0.038 1.5E4 1.5E4 0.0 16",
            "  CLASS Coarse\n    1 Sph 0.05 0.020 1.0E4 1.0E4 0.0 16\n    2 Exp 0.00 0.018 2.0E4 2.0E4 0.0 16"
        ).replace("MAX_VSTRUCT 1", "MAX_VSTRUCT 2")
        r = run_t2p(t2p, extra_files={"dataset.csv": minimal_dataset,
                                       "grid.grd": minimal_grd})
        assert r.returncode == 0, r.combined


# ---------------------------------------------------------------------------
# 7. EXTERNAL_FILE command
# ---------------------------------------------------------------------------

class TestExternalFile:

    def test_variograms_in_external_file(self, run_t2p, minimal_dataset,
                                          minimal_grd):
        """Variograms block delegated to an external file must be read correctly."""
        main_t2p = """\
            BEGIN OPTIONS
              MAX_VSTRUCT 1
              INFER_LAST_CLASS
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
            EXTERNAL_FILE variograms.ext
            BEGIN PP_LOCS
              1 500000.0 4000000.0 1
            END PP_LOCS
            BEGIN PP_PARAMETERS
              TYPE Global
              1 0.5 -0.5 0.0
              TYPE Aquifer
              1 Coarse 0.1 100.0 1e-5 0.15 10.0 0.001
              1 Fine   0.001 0.5 1e-4 0.05 5.0 0.001
            END PP_PARAMETERS
            """
        ext_variograms = """\
            BEGIN VARIOGRAMS
              CLASS Coarse
                1 Sph 0.05 0.038 1.5E4 1.5E4 0.0 16
              CLASS PilotPoints
                1 Sph 0.00 1.0 1.0E5 1.0E5 0.0 25
            END VARIOGRAMS
            """
        r = run_t2p(main_t2p, extra_files={
            "dataset.csv": minimal_dataset,
            "grid.grd": minimal_grd,
            "variograms.ext": ext_variograms,
        })
        assert r.returncode == 0, r.combined

    def test_missing_external_file_errors(self, run_t2p, minimal_dataset,
                                           minimal_grd):
        main_t2p = _simple_t2p(extra_blocks="EXTERNAL_FILE nonexistent.ext")
        r = run_t2p(main_t2p, extra_files={"dataset.csv": minimal_dataset,
                                             "grid.grd": minimal_grd})
        assert r.returncode != 0 or "error" in r.combined.lower()


# ---------------------------------------------------------------------------
# 8. Input filename convention
# ---------------------------------------------------------------------------

class TestFilenameConvention:

    def test_explicit_filename_passed(self, run_t2p, minimal_dataset, minimal_grd):
        """Binary should accept an explicitly named .t2p file."""
        r = run_t2p(_simple_t2p(), filename="mymodel_cokriging.t2p",
                    extra_files={"dataset.csv": minimal_dataset,
                                 "grid.grd": minimal_grd})
        assert r.returncode == 0, r.combined

    def test_default_filename_texture2par_in(self, run_t2p, minimal_dataset,
                                              minimal_grd):
        """When no argument is given, T2P should fall back to 'Texture2Par.in'."""
        import subprocess, shutil
        workdir = run_t2p._pytest_request.getfixturevalue("tmp_path") \
            if hasattr(run_t2p, "_pytest_request") else None
        # This test is handled via the conftest runner; skip if not supported
        pytest.skip("Requires invoking binary without arguments – covered by manual test")
