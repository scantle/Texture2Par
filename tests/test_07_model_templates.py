"""Regression coverage for template-driven groundwater-model output files."""

from __future__ import annotations

from .test_02_outputs import _base_t2p
from .test_03_aquitard_pilot_points import COMMON_FILES as IWFM_FILES
from .test_03_aquitard_pilot_points import _input as _iwfm_input
from .test_03_interpolation import (
    BLOCK_FILES,
    _block_modflow_t2p,
)


MODFLOW_LPF_TEMPLATE = """\
0 -1.0E+30 0
0
0
1.0
0
0
CONSTANT 1.0
CONSTANT 1.0
"""


IWFM_GROUNDWATER_TEMPLATE = """\
4.0
boundary.dat
pumping.dat
subsidence.dat
# Debugging
0
# Hydrograph output count
0
# Element-face output count
0
# Aquifer parameters
0
# Conversion factors
1.0
# Units
unit
# Node parameters
0
0
0
0
0
0
"""


IWFM_SUBSIDENCE_TEMPLATE = """\
V4
# Input and output files
0
# Subsidence output data
0
# Print settings
0
# Subsidence parameters
0
# Conversion factors
1.0 1.0 1.0 1.0 0.1 0.2
# Node parameters
1 0.0 0.0 0.0 0.1 0.2
2 0.0 0.0 0.0 0.1 0.2
3 0.0 0.0 0.0 0.1 0.2
4 0.0 0.0 0.0 0.1 0.2
5 0.0 0.0 0.0 0.1 0.2
6 0.0 0.0 0.0 0.1 0.2
"""


def _modflow_parameter_input(template="model_template.lpf") -> str:
    text = _base_t2p().replace(
        "TYPE GRID\n          GRID_FILE grid.grd",
        "TYPE MODFLOW\n"
        "          NAM_FILE model.nam\n"
        f"          TEMPLATE_FILE {template}",
    )
    return text.replace("5050.0 4005.0", "1.0 1.0")


def _iwfm_parameter_input(template="groundwater_template.dat") -> str:
    text = _iwfm_input(explicit_zones=False).replace("  WRITE_NONE\n", "")
    return text.replace(
        "  PREPROC_FILE preprocessor.in",
        "  PREPROC_FILE preprocessor.in\n"
        f"  TEMPLATE_FILE {template}",
    )


def _iwfm_subsidence_input(template="subsidence_template.dat") -> str:
    text = _iwfm_parameter_input().replace(
        "BEGIN PP_LOCS",
        "BEGIN SUBSIDENCE\n"
        "  FINE_CLASSES Fine\n"
        f"  TEMPLATE_FILE {template}\n"
        "  VARIOGRAM\n"
        "    1 Sph 0.0 1.0 100.0 100.0 0.0 1\n"
        "END SUBSIDENCE\n\n"
        "BEGIN PP_LOCS",
    )
    return text.replace(
        "100.0 100.0 1e-5 0.15 1.0 0.0",
        "100.0 100.0 0.01 1e-6 0.0 0.15 1.0 0.0",
    ).replace(
        "0.01 0.01 1.0 0.0",
        "0.01 0.01 0.001 1e-5 1e-4 1.0 0.0",
    )


class TestModflowTemplateOutput:

    def test_missing_lpf_target_is_created_from_template(self, run_t2p):
        files = {
            **{name: value for name, value in BLOCK_FILES.items()
               if name != "model.lpf"},
            "model_template.lpf": MODFLOW_LPF_TEMPLATE,
        }
        result = run_t2p(_modflow_parameter_input(), extra_files=files)
        assert result.returncode == 0, result.combined
        output = result.output_file("model.lpf")
        assert output.exists()
        assert "Layer Property Flow (LPF) Package Written by Texture2Par" in (
            output.read_text(encoding="utf-8")
        )

    def test_missing_modflow_template_has_clear_error(self, run_t2p):
        files = {name: value for name, value in BLOCK_FILES.items()
                 if name != "model.lpf"}
        result = run_t2p(
            _modflow_parameter_input("missing_template.lpf"),
            extra_files=files,
        )
        assert result.returncode != 0
        assert "MODFLOW LPF/UPW template file does not exist" in result.combined
        assert "missing_template.lpf" in result.combined

    def test_only_texture_needs_neither_target_nor_template(self, run_t2p):
        files = {name: value for name, value in BLOCK_FILES.items()
                 if name != "model.lpf"}
        result = run_t2p(_block_modflow_t2p(), extra_files=files)
        assert result.returncode == 0, result.combined
        assert not result.output_exists("model.lpf")


class TestIwfmTemplateOutput:

    def test_missing_groundwater_target_is_created_from_template(self, run_t2p):
        result = run_t2p(
            _iwfm_parameter_input(),
            extra_files={
                **IWFM_FILES,
                "groundwater_template.dat": IWFM_GROUNDWATER_TEMPLATE,
            },
        )
        assert result.returncode == 0, result.combined
        output = result.output_file("groundwater.dat")
        assert output.exists()
        assert output.stat().st_size > 0

    def test_missing_groundwater_template_has_clear_error(self, run_t2p):
        result = run_t2p(
            _iwfm_parameter_input("missing_groundwater_template.dat"),
            extra_files=IWFM_FILES,
        )
        assert result.returncode != 0
        assert "IWFM groundwater template file does not exist" in result.combined
        assert "missing_groundwater_template.dat" in result.combined

    def test_write_none_needs_no_groundwater_template(self, run_t2p):
        result = run_t2p(
            _iwfm_input(explicit_zones=False),
            extra_files=IWFM_FILES,
        )
        assert result.returncode == 0, result.combined
        assert not result.output_exists("groundwater.dat")

    def test_missing_subsidence_target_is_created_from_template(self, run_t2p):
        result = run_t2p(
            _iwfm_subsidence_input(),
            extra_files={
                **IWFM_FILES,
                "groundwater_template.dat": IWFM_GROUNDWATER_TEMPLATE,
                "subsidence_template.dat": IWFM_SUBSIDENCE_TEMPLATE,
            },
        )
        assert result.returncode == 0, result.combined
        output = result.output_file("subsidence.dat")
        assert output.exists()
        assert output.stat().st_size > 0

    def test_missing_subsidence_template_has_clear_error(self, run_t2p):
        result = run_t2p(
            _iwfm_subsidence_input("missing_subsidence_template.dat"),
            extra_files={
                **IWFM_FILES,
                "groundwater_template.dat": IWFM_GROUNDWATER_TEMPLATE,
            },
        )
        assert result.returncode != 0
        assert "IWFM subsidence template file does not exist" in result.combined
        assert "missing_subsidence_template.dat" in result.combined
