"""
conftest.py – shared fixtures for the Texture2Par v2 test suite.

Place Texture2Par.exe (or a shell wrapper named texture2par) on PATH,
or set the T2P_EXE environment variable to the full path before running.

    export T2P_EXE=/path/to/Texture2Par.exe
    pytest tests/ -v
"""

import os
import shutil
import subprocess
import textwrap
from pathlib import Path
from .data import aem2d, grid2d, gridspec, pc2d
import pytest


# ---------------------------------------------------------------------------
# Locate the binary
# ---------------------------------------------------------------------------

def _find_exe() -> str:
    if os.path.exists("bin/Texture2Par.exe"):
        return os.path.abspath("bin/Texture2Par.exe")
    exe = os.environ.get("T2P_EXE")
    if exe:
        return exe
    for name in ("Texture2Par.exe", "Texture2Par", "texture2par"):
        found = shutil.which(name)
        if found:
            return found
    pytest.skip("Texture2Par binary not found. Set T2P_EXE env var.")


T2P_EXE = _find_exe()


# ---------------------------------------------------------------------------
# Core runner fixture
# ---------------------------------------------------------------------------

@pytest.fixture
def run_t2p(tmp_path):
    """
    Returns a callable: run_t2p(input_file, extra_files=None)
    Copies every file in extra_files to tmp_path, then runs Texture2Par
    with input_file (also written to tmp_path).

    Returns a CompletedProcess-like namespace:
        .returncode  int
        .stdout      str
        .stderr      str
        .log         str   (contents of <input_stem>.out, if present)
        .workdir     Path
    """
    def _run(t2p_text: str, filename: str = "test.t2p",
             extra_files: dict[str, str] | None = None) -> "RunResult":
        workdir = tmp_path
        (workdir / filename).write_text(textwrap.dedent(t2p_text))
        if extra_files:
            for fname, content in extra_files.items():
                (workdir / fname).write_text(textwrap.dedent(content))

        result = subprocess.run(
            [T2P_EXE, filename],
            cwd=workdir,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
        )

        log_path = workdir / (Path(filename).stem + ".out")
        log = log_path.read_text(encoding="utf-8", errors="replace") if log_path.exists() else ""

        return RunResult(
            returncode=result.returncode,
            stdout=result.stdout,
            stderr=result.stderr,
            log=log,
            workdir=workdir,
        )

    return _run


class RunResult:
    def __init__(self, returncode, stdout, stderr, log, workdir):
        self.returncode = returncode
        self.stdout = stdout
        self.stderr = stderr
        self.log = log
        self.workdir = workdir

    @property
    def combined(self):
        return self.stdout + self.stderr + self.log

    def output_file(self, name: str) -> Path:
        return self.workdir / name

    def output_exists(self, name: str) -> bool:
        return (self.workdir / name).exists()


# ---------------------------------------------------------------------------
# Minimal valid dataset fixture (reused across many tests)
# ---------------------------------------------------------------------------

MINIMAL_DATASET = """\
Location,WellIndex,IntervalIndex,X,Y,Zland,TopDepth,BotDepth,Coarse,Fine
W1,1,1,5000.0,4000.0,100.0,0.0,5.0,0.6,0.4
W1,1,2,5000.0,4000.0,100.0,5.0,10.0,0.5,0.5
W2,2,1,5100.0,4000.0,100.0,0.0,5.0,0.7,0.3
W2,2,2,5100.0,4000.0,100.0,5.0,10.0,0.4,0.6
W3,3,1,5050.0,4010.0,100.0,0.0,5.0,0.55,0.45
W3,3,2,5050.0,4010.0,100.0,5.0,10.0,0.45,0.55
"""

MINIMAL_GRD = """\
NX    5
NY    5
NZ    2
DX    2000.0
DY    2000.0
TOP   constant 100.0
BOT   constant 80.0
IBOUND constant 1
"""


@pytest.fixture
def minimal_dataset():
    return MINIMAL_DATASET


@pytest.fixture
def minimal_grd():
    return MINIMAL_GRD


@pytest.fixture
def butte_ok_data():
    return {
        "kws":dict(
            datafile="pc2d.csv",
            gridfile="grid.grd",
            classes="Coarse\n  Fine",
            variogram="sph 0.0 0.12 5000 5000 0.0 100",
        ),
        "files":{
            "pc2d.csv":pc2d.pc2d,
            "grid.grd":gridspec.gridspec,
        }
    }

@pytest.fixture
def butte_ck_data():
    return {
        "kws":dict(
            datafile=["pc2d.csv", "aem2d.csv"],
            gridfile="grid.grd",
            classes="Coarse\n  Fine",
            covariate="Rho",
            variogram="sph 0.0 0.12 5000 5000 0.0 100",
        ),
        "files":{
            "pc2d.csv":pc2d.pc2d,
            "aem2d.csv":aem2d.aem2d,
            "grid.grd":gridspec.gridspec,
        }
    }