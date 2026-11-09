"""Focused tests for the contained resistivity-to-texture preprocessor."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest


GRID = """\
NX    2
NY    2
NZ    1
DX    1.0
DY    1.0
TOP   constant 100.0
DZ    10.0
IBOUND constant 1
"""


DATASET = """\
Location,WellIndex,IntervalIndex,X,Y,Zland,TopDepth,BotDepth,LogRho,var_LogRho
W1,1,1,0.1,0.1,100.0,0.0,10.0,3.5,0.0
W2,2,1,1.9,0.1,100.0,0.0,10.0,4.0,0.0
W3,3,1,0.1,1.9,100.0,0.0,10.0,4.5,0.0
W4,4,1,1.9,1.9,100.0,0.0,10.0,5.0,0.0
"""


SCALES = """\
node,layer,Fine_scale,Coarse_scale,Fine_shape,Coarse_shape
1,1,20.0,100.0,0.25,0.25
2,1,20.0,100.0,0.25,0.25
3,1,20.0,100.0,0.25,0.25
4,1,20.0,100.0,0.25,0.25
"""


def _multipliers(value: float) -> str:
    return "node,layer,multiplier\n" + "".join(
        f"{node},1,{value}\n" for node in range(1, 5)
    )


def _preproc_t2p(mult_file: str = "mult.csv") -> str:
    return f"""\
        BEGIN OPTIONS
          MAX_VSTRUCT 1
          INTERP_METHOD OK
          ONLY_TEXTURE
          WRITE_NODE_FILES
          WRITE_VARIANCE
          KVME
        END OPTIONS

        BEGIN FLOW_MODEL
          TYPE GRID
          GRID_FILE grid.grd
        END FLOW_MODEL

        BEGIN CLASSES
          Fine
          Coarse
        END CLASSES

        BEGIN PREPROC
          MODE RES2TEXTURE
          CLASS LogRho
          SCALE_FILE scales.csv
          MULT_FILE {mult_file}
        END PREPROC

        BEGIN DATASET
          FILE dataset.csv
        END DATASET

        BEGIN VARIOGRAMS
          CLASS LogRho
            1 Sph 0.0 1.0 10.0 10.0 0.0 4
        END VARIOGRAMS
        """


def _output(workdir: Path, project: str, name: str) -> pd.DataFrame:
    path = workdir / f"{project}_{name}.csv"
    assert path.exists(), sorted(item.name for item in workdir.iterdir())
    return pd.read_csv(path)


def _files(multiplier: float, scales: str = SCALES) -> dict[str, str]:
    return {
        "grid.grd": GRID,
        "dataset.csv": DATASET,
        "scales.csv": scales,
        "mult.csv": _multipliers(multiplier),
    }


def test_res2texture_probabilities_are_normalized(run_t2p):
    result = run_t2p(_preproc_t2p(), filename="normalized.t2p",
                     extra_files=_files(0.0))

    assert result.returncode == 0, result.combined
    fine = _output(result.workdir, "normalized", "FINE")["Layer1"].to_numpy()
    coarse = _output(result.workdir, "normalized", "COARSE")["Layer1"].to_numpy()
    assert np.isfinite(fine).all()
    assert np.isfinite(coarse).all()
    assert fine + coarse == pytest.approx(np.ones(4), abs=2e-6)


def test_res2texture_multiplier_uses_kriging_uncertainty(run_t2p):
    zero = run_t2p(_preproc_t2p(), filename="zero.t2p",
                   extra_files=_files(0.0))
    assert zero.returncode == 0, zero.combined
    zero_fine = _output(zero.workdir, "zero", "FINE")["Layer1"].to_numpy().copy()
    zero_raw = _output(zero.workdir, "zero", "LOGRHO").copy()

    shifted = run_t2p(_preproc_t2p(), filename="shifted.t2p",
                      extra_files=_files(1.0))
    assert shifted.returncode == 0, shifted.combined
    shifted_fine = _output(shifted.workdir, "shifted", "FINE")["Layer1"].to_numpy()
    shifted_raw = _output(shifted.workdir, "shifted", "LOGRHO")

    assert (zero_raw["Var_Layer1"] > 0.0).all()
    assert shifted_raw["Layer1"].to_numpy() == pytest.approx(
        zero_raw["Layer1"].to_numpy(), abs=1e-10
    )
    assert shifted_raw["Var_Layer1"].to_numpy() == pytest.approx(
        zero_raw["Var_Layer1"].to_numpy(), abs=1e-10
    )
    assert not np.allclose(shifted_fine, zero_fine, atol=1e-8, rtol=0.0)


def test_res2texture_rejects_incomplete_scale_grid(run_t2p):
    incomplete = "\n".join(SCALES.splitlines()[:-1]) + "\n"
    result = run_t2p(_preproc_t2p(), filename="incomplete.t2p",
                     extra_files=_files(0.0, scales=incomplete))

    assert result.returncode != 0
    assert "does not cover every model cell" in result.combined
