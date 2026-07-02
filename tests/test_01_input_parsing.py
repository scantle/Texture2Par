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

    def test_pp_ib_interp_method_defaults_are_ok(
        self, run_t2p, minimal_dataset, minimal_grd
    ):
        r = run_t2p(
            _simple_t2p(),
            extra_files={"dataset.csv": minimal_dataset, "grid.grd": minimal_grd},
        )
        assert r.returncode == 0, r.combined
        pp_line = next(
            line for line in r.log.splitlines()
            if "Pilot-Point Interpolation Method:" in line
        )
        ib_line = next(
            line for line in r.log.splitlines()
            if "Interbed Interpolation Method:" in line
        )
        assert "Ordinary Kriging (libkriging)" in pp_line
        assert "Ordinary Kriging (libkriging)" in ib_line

    @pytest.mark.parametrize("option", ["PP_INTERP_METHOD SK", "IB_INTERP_METHOD SK"])
    def test_pp_ib_interp_method_sk_is_rejected(
        self, run_t2p, minimal_dataset, minimal_grd, option
    ):
        r = run_t2p(
            _simple_t2p(options_extra=option),
            extra_files={"dataset.csv": minimal_dataset, "grid.grd": minimal_grd},
        )
        assert r.returncode != 0
        assert option.split()[0] in r.combined
        assert "defaults to OK" in r.combined
        assert "SK mean handling is not defined" in r.combined

    @pytest.mark.parametrize("option", [
        "WRITE_DATASET_FILES",
        "WRITE_NODE_FILES",
        "WRITE_ELEVATION_FILE",
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

    @pytest.mark.parametrize("mode", ["LIBKRIGING", "LEGACY_STANDARDIZED"])
    def test_ck_mean_mode_accepted_and_reported(
        self, run_t2p, minimal_dataset, minimal_grd, mode
    ):
        r = run_t2p(_simple_t2p(options_extra=f"CK_MEAN_MODE {mode}"),
                    extra_files={"dataset.csv": minimal_dataset,
                                 "grid.grd": minimal_grd})
        assert r.returncode == 0, r.combined
        summary_line = next(
            line for line in r.log.splitlines() if "CK Mean Mode:" in line
        )
        assert summary_line.rstrip().endswith(mode)

    def test_unknown_ck_mean_mode_rejected(
        self, run_t2p, minimal_dataset, minimal_grd
    ):
        r = run_t2p(_simple_t2p(options_extra="CK_MEAN_MODE UNKNOWN"),
                    extra_files={"dataset.csv": minimal_dataset,
                                 "grid.grd": minimal_grd})
        assert r.returncode != 0
        assert "Invalid CK_MEAN_MODE: UNKNOWN" in r.combined


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

    def test_missing_grid_file_errors(self, run_t2p, minimal_dataset):
        """GRID_FILE referencing a non-existent file must error."""
        r = run_t2p(_simple_t2p(),
                    extra_files={"dataset.csv": minimal_dataset})
        assert r.returncode != 0 or "error" in r.combined.lower()


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
