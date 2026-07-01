# Texture2Par v2 – Automated Test Suite

A pytest-based test suite that exercises the Texture2Par v2 binary end-to-end.
Each test writes real `.t2p` input files to a temporary directory, runs the
binary, and asserts on its return code, log output, and output files.

---

## Structure

```
tests/
  conftest.py               # shared fixtures (binary locator, run_t2p runner)
  test_01_input_parsing.py  # input syntax, required fields, and validation
  test_02_outputs.py        # output files, columns, elevation, and subsidence
  test_03_interpolation.py  # core interpolation, 3D OK, SGS, and weight caching
  test_04_butte_model.py    # larger Butte ordinary-kriging reference tests
  test_05_butte_model_ck.py # larger Butte cokriging and 3D CK reference tests
pytest.ini               # pytest configuration
requirements-test.txt    # pytest dependency
```

The numeric prefixes make the default pytest run grow from basic parsing and
output contracts through interpolation behavior to the larger Butte reference
cases. They control presentation order only; tests remain independent.

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

`T2P_EXE` takes precedence over both the repository's `bin/Texture2Par.exe`
and executables found on `PATH`. Pytest prints the resolved executable path in
its session header so regression baselines can be audited.

**Option B – add the binary to PATH**  
The suite searches for `Texture2Par.exe`, `Texture2Par`, and `texture2par` in
that order.

---

## Running

```bash
# Run everything
pytest

# If pytest cannot access its default temp directory on Windows
pytest --basetemp=.pytest_tmp

# Run a single file
pytest tests/test_01_input_parsing.py

# Run a single class or test
pytest tests/test_03_interpolation.py::TestSGS
pytest tests/test_02_outputs.py::TestSubsidenceOutputs::test_ib_thk_non_negative

# Stop on first failure
pytest -x

# Show full output on failure
pytest --tb=long
```

---

## Test categories

### `test_01_input_parsing.py`

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

### `test_02_outputs.py`

| Class | What it checks |
|---|---|
| `TestLogFile` | log created; elapsed time; options echoed; errors on bad run |
| `TestNodeFiles` | KH/KV/SS/SY written; spatial columns present; KH > 0; SY in [0,1]; suppressed without flag |
| `TestDatasetFiles` | Coarse file written; X/Y columns present; suppressed without flag |
| `TestElevationFile` | written with flag; suppressed without flag |
| `TestWriteNone` | model output files suppressed |
| `TestSubsidenceOutputs` | run completes; IB_THK/KVF written; IB_THK ≥ 0 |
| `TestOutputRowCounts` | KH row count matches grid size (5×5×2 = 50) |

### `test_03_interpolation.py`

| Class | What it checks |
|---|---|
| `TestOrdinaryKriging` | run completes; KH written; texture bounded [0,1]; INFER classes |
| `TestLibKrigingOrdinary` | optional OK_LIB data and pilot-point interpolation; 2D/3D isotropic parity; explicit libkriging/GSLIB anisotropy and anisotropic-search coverage; NNEAR; nested structures; supported variogram mappings; duplicate and unsupported-feature diagnostics |
| `TestCorrectWeights` | CORRECT_WEIGHTS flag |
| `TestCoincidentLocations` | duplicate XY/XYZ conditioning points are fatal; exact target hits remain valid |
| `TestSGS` | flag accepted; different seeds differ; same seed reproducible; SGS + weight write incompatibility |
| `TestWeightCaching` | WRITE/READ_GRID_WEIGHTS; WRITE/READ_PP_WEIGHTS; TEXT suffix |
| `TestOnlyTexture` | no pilot points required; no hydraulic property files written |
| `TestMaxLogLength` | INTERP_DIM 3 log interval splitting and processed XYZ support points |
| `TestInterpDim3GridGeometry` | 3D grid elevation output and implied layer-center Z coordinates |
| `TestOrdinaryKriging3D` | 3D OK exact hits; multi-log reference checks; rotation/anisotropy; nnear; ANISOTROPIC_SEARCH; Exp/Gau |
| `TestReadGridded` | pre-interpolated grid input accepted |
| `TestAnisotropicSearch` | ANISOTROPIC_SEARCH flag |
| `TestNoTexBounds` | NO_TEX_BOUNDS allows values outside [0,1] |

### `test_04_butte_model.py` and `test_05_butte_model_ck.py`

The Butte tests provide larger ordinary-kriging and cokriging reference cases,
including anisotropy, 3D interpolation, and CK covariance checks. They run last
because they are broader and more computationally involved than the basic
contract tests above.

### Cokriging Reference Tests

Butte CK reference tests use a full augmented CK matrix in Python. For
secondary data, that Python reference uses the entire selected secondary
subsample for every target. Texture2Par's default `MAX_SECONDARY_DIST` is
finite, so tests that compare to the full Python reference must explicitly set
a very large `MAX_SECONDARY_DIST`, such as `1e10`. Tests for default or
restricted secondary-distance behavior should use their own reference logic or
be treated as separate behavior/regression tests.

The 3D OK and CK Python references intentionally follow Texture2Par's **legacy**
3D anisotropy convention, including the positive-downward dip / Fortran
rotation-matrix convention. This convention is unconventional and differs from
libkriging/scipy/GSLIB-style rotation. The current expected values, including
those for `ANISOTROPIC_SEARCH`, must not be updated casually. They are expected
to change deliberately during the libkriging migration, with the convention
change reviewed as such. A compact traced-pair covariance regression guards
the CK reference against accidental dip-sign changes in the meantime.

`OK_LIB` parity tests are therefore deliberately isotropic. `OK_LIB` passes
angles to libkriging using its GSLIB-style rotation convention, while legacy
`OK` retains Texture2Par's historical convention. Any anisotropic `OK_LIB`
test must state which convention it validates and must not reuse legacy
anisotropic expectations without explicitly documenting the expected
difference. The two sets of anisotropic tests must not be merged into one
parity expectation. Cutting over default `OK` to libkriging will intentionally
change anisotropic results where the historical rotation convention differs.

`ANISOTROPIC_SEARCH` controls only the `OK_LIB` neighbor-selection metric;
covariance evaluation continues to use every structure's own variogram
anisotropy regardless of that flag. Search currently uses the first variogram
structure's major/minor range ratios and angles. Sector search remains a TODO
and is not exposed by Texture2Par in this migration stage.

At this migration stage, `OK_LIB` is optional and applies to dataset,
pilot-point, and interbed ordinary kriging. The pilot-point and interbed paths
pack geometry once per category and interpolate multiple value columns through
the libkriging adapter. The first solve stores libkriging's in-memory weights;
later columns update only the observation values and reuse those weights.
Legacy PP/IB weight files are intentionally unsupported under `OK_LIB`; normal
`OK` continues to use its existing dense-weight paths.

In `INTERP_DIM 3`, processed observations remain XYZ support points; they are
not reduced to the per-model-layer XY averages used by 2D layer interpolation.
Grid target Z coordinates are model-layer center elevations. Duplicate
conditioning support coordinates are fatal in either dimension, while a grid
target coincident with one observation is a valid exact hit.

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
