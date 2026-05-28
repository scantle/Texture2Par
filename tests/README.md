# Texture2Par v2 – Automated Test Suite

A pytest-based test suite that exercises the Texture2Par v2 binary end-to-end.
Each test writes real `.t2p` input files to a temporary directory, runs the
binary, and asserts on its return code, log output, and output files.

---

## Structure

```
tests/
  conftest.py            # shared fixtures (binary locator, run_t2p runner)
  test_input_parsing.py  # OPTIONS / CLASSES / DATASET / FLOW_MODEL / VARIOGRAMS blocks
  test_interpolation.py  # ordinary kriging, SGS, weight caching, ONLY_TEXTURE, READ_GRIDDED
  test_outputs.py        # log file, node/dataset CSVs, elevation file, subsidence outputs
pytest.ini               # pytest configuration
requirements-test.txt    # pytest dependency
```

---

## Setup

```bash
pip install -r requirements-test.txt
```

Point the suite at the binary in one of two ways:

**Option A – environment variable (recommended)**
```bash
export T2P_EXE=/path/to/Texture2Par.exe   # Linux / macOS
set T2P_EXE=C:\path\to\Texture2Par.exe    # Windows
```

**Option B – add the binary to PATH**  
The suite searches for `Texture2Par.exe`, `Texture2Par`, and `texture2par` in
that order.

---

## Running

```bash
# Run everything
pytest

# Run a single file
pytest tests/test_input_parsing.py

# Run a single class or test
pytest tests/test_interpolation.py::TestSGS
pytest tests/test_outputs.py::TestSubsidenceOutputs::test_ib_thk_non_negative

# Stop on first failure
pytest -x

# Show full output on failure
pytest --tb=long
```

---

## Test categories

### `test_input_parsing.py`

| Class | What it checks |
|---|---|
| `TestOptionsBlock` | MAX_VSTRUCT mandatory; boolean flags; INTERP_DIM; NO_DATA_VALUE; PROJECT_NAME; SEED |
| `TestClassesBlock` | single class + INFER; underscore/hyphen names; rejection of invalid characters |
| `TestDatasetBlock` | missing file detection; all 7 TRANSFORM keywords; -999 no-data rows |
| `TestFlowModelBlock` | GRID type acceptance; unknown TYPE rejection; missing GRID_FILE |
| `TestGrdFile` | constant TOP/BOT; DZ instead of BOT; USE_UNIFORM_IBOUND; missing NX |
| `TestVariogramsBlock` | all 5 variogram types; PilotPoints must be last; missing PP variogram; composite structures |
| `TestExternalFile` | variogram block in external file; missing external file |
| `TestFilenameConvention` | explicit .t2p filename |

### `test_interpolation.py`

| Class | What it checks |
|---|---|
| `TestOrdinaryKriging` | run completes; KH written; texture bounded [0,1]; INFER classes |
| `TestCorrectWeights` | CORRECT_WEIGHTS flag |
| `TestSGS` | flag accepted; different seeds differ; same seed reproducible; SGS + weight write incompatibility |
| `TestWeightCaching` | WRITE/READ_GRID_WEIGHTS; WRITE/READ_PP_WEIGHTS; TEXT suffix |
| `TestOnlyTexture` | no pilot points required; no hydraulic property files written |
| `TestReadGridded` | pre-interpolated grid input accepted |
| `TestAnisotropicSearch` | ANISOTROPIC_SEARCH flag |
| `TestNoTexBounds` | NO_TEX_BOUNDS allows values outside [0,1] |

### `test_outputs.py`

| Class | What it checks |
|---|---|
| `TestLogFile` | log created; elapsed time; options echoed; errors on bad run |
| `TestNodeFiles` | KH/KV/SS/SY written; spatial columns present; KH > 0; SY in [0,1]; suppressed without flag |
| `TestDatasetFiles` | Coarse file written; X/Y columns present; suppressed without flag |
| `TestElevationFile` | written with flag; suppressed without flag |
| `TestWriteNone` | model output files suppressed |
| `TestSubsidenceOutputs` | run completes; IB_THK/KVF written; IB_THK ≥ 0 |
| `TestOutputRowCounts` | KH row count matches grid size (5×5×2 = 50) |

---

## Adding new tests

All tests follow the same pattern:

```python
def test_my_new_feature(self, run_t2p):
    r = run_t2p(
        "<your .t2p content here>",
        extra_files={"dataset.csv": DATASET, "grid.grd": GRD},
    )
    assert r.returncode == 0, r.combined   # r.combined = stdout + stderr + log
    assert r.output_exists("some_output.csv")
```

`run_t2p` writes all files to a fresh `tmp_path`, runs the binary, and returns
a `RunResult` with `.returncode`, `.stdout`, `.stderr`, `.log`, `.combined`,
`.workdir`, `.output_file(name)`, and `.output_exists(name)`.
