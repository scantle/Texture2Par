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
| `TestLibKrigingOrdinary` | default OK data, pilot-point, and interbed routing; OK_LIB alias; rejection of removed OK_LEGACY; explicit libkriging/GSLIB anisotropy and anisotropic-search coverage; NNEAR; nested structures; supported variogram mappings; duplicate and unsupported-feature diagnostics |
| `TestLibKrigingSimple` | default SK routing through the shared libkriging backend; global-mean 2D/3D references; SK/OK distinction; retained legacy-OK pilot-point route; unsupported-feature diagnostics |
| `TestCorrectWeights` | CORRECT_WEIGHTS flag |
| `TestCoincidentLocations` | duplicate XY/XYZ conditioning points are fatal; exact target hits remain valid |
| `TestSGS` | flag accepted; different seeds differ; same seed reproducible; SGS + weight write incompatibility |
| `TestWeightCaching` | clear unsupported diagnostics for external grid and pilot-point saved weights |
| `TestOnlyTexture` | no pilot points required; no hydraulic property files written |
| `TestMaxLogLength` | INTERP_DIM 3 log interval splitting and processed XYZ support points |
| `TestInterpDim3GridGeometry` | 3D grid elevation output and implied layer-center Z coordinates |
| `TestOrdinaryKriging3D` | libkriging 3D OK exact hits and reference checks; rotation/anisotropy; nnear; ANISOTROPIC_SEARCH; Exp/Gau |
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

Historical 3D helper tests retain explicit coverage of Texture2Par's old
positive-downward dip / Fortran rotation matrix, but candidate executable
interpolation tests follow libkriging's GSLIB-style convention. Historical
candidate-output parity is no longer routed through a legacy deterministic
backend.

Default `INTERP_METHOD OK` now uses libkriging for ordinary kriging and passes
angles using libkriging's GSLIB-style rotation convention. `OK_LIB` remains a
temporary alias for that same path. `OK_LEGACY` has been removed. Any
anisotropic default-OK test must state which convention it
validates and must not reuse legacy
anisotropic expectations without explicitly documenting the expected
difference. The two sets of anisotropic tests must not be merged into one
parity expectation. The default-OK cutover intentionally changes anisotropic
results where the historical rotation convention differs.

`ANISOTROPIC_SEARCH` controls only the libkriging OK/SK neighbor-selection metric;
covariance evaluation continues to use every structure's own variogram
anisotropy regardless of that flag. Search currently uses the first variogram
structure's major/minor range ratios and angles. Sector search remains a TODO
and is not exposed by Texture2Par in this migration stage.

`INTERP_METHOD` controls model/data interpolation. Pilot-point and interbed
interpolation are configured independently with `PP_INTERP_METHOD` and
`IB_INTERP_METHOD`; both default to `OK` and also accept `IDW`. Their `OK`
setting uses libkriging, regardless of whether model/data interpolation uses
`OK`, `OK_LIB`, `SK`, `IDW`, or `OK + SGSIM`. `SK` is rejected for PP/interbed
because those workflows do not define an SK mean convention.

The pilot-point and interbed libkriging paths pack geometry
once per category and interpolate multiple value columns through the
libkriging adapter. The first solve stores libkriging's in-memory weights;
later columns update only the observation values and reuse those weights.
External grid, pilot-point, and interbed saved-weight files remain unsupported
under their respective libkriging OK paths and fail with a diagnostic. `OK + SGSIM` continues to select the existing
SGSIM implementation; `OK_LIB + SGSIM` is rejected because the alias names the
ordinary-kriging backend directly.

Default `INTERP_METHOD SK` now uses the same `t_libkrige` backend as OK. The
backend mode selects `unbias=1` for OK and `unbias=0` for SK; SK passes the
dataset/class global mean to libkriging's `set_obs`. That mean is computed by
the existing dataset readers across all layers/categories (and is
thickness-weighted for 3D interval data), not recomputed per category or
target. SK uses the libkriging search path and `NNEAR`, whereas legacy
the removed `t_globkrige` assembled all observations in the current category.

SK follows libkriging's GSLIB-style anisotropy convention, not Texture2Par's
historical rotation behavior. Because PP/interbed methods default independently
to OK, `INTERP_METHOD SK` no longer routes either workflow through legacy
ordinary kriging. External saved weights, cokriging, SGSIM, block kriging, and
`INTERP_ORDER AFTER` are rejected for SK rather than silently falling back;
the latter generated parameter data never had an initialized legacy SK mean.

Deterministic secondary-data interpolation under default `OK` and its `OK_LIB`
alias now uses the same libkriging backend with `nvar > 1`. Observations are
packed separately for the primary and each connected secondary variable.
Both CK mean modes use the same single-constraint libkriging covariance,
search, and solve path; only value preparation and libkriging's existing mean
correction flag differ. Only the primary estimate and variance are scattered
back to the model.

`CK_MEAN_MODE LIBKRIGING` is the default. Secondary values are passed without
Texture2Par's historical global shift, and libkriging applies its non-standard
CK target-neighbourhood mean correction. `CK_MEAN_MODE LEGACY_STANDARDIZED`
instead shifts every secondary by `primary_global_mean-secondary_global_mean`
and disables the local correction. This preserves Texture2Par's historical
standardized mean formulation while still using libkriging. The underlying libkriging
`nonstd_ck_mean_correction` option remains true by default for other callers.

Each auto- and cross-variogram in the zoned `global_vario_idx` map is supplied
to libkriging for every variable pair, including all configured nested
structures. `NNEAR` and search anisotropy are applied per variable from that
variable's auto-variogram. `MAX_PRIMARY_DIST` applies to the primary variable;
`MAX_SECONDARY_DIST` applies to every secondary variable. When the latter is
zero, the retained Texture2Par default is three times the largest major range
in the secondary auto-variogram. Full-secondary Python references must still
set a large explicit value such as `1e10`.

New-path CK follows libkriging/GSLIB-style covariance and search anisotropy.
SK cokriging, asymmetric cross-variograms, block
kriging, sector search, SGSIM through the direct libkriging alias, and external
saved-weight files remain unsupported. IDW and the existing `OK + SGSIM`
routing are unchanged. `kriging.f90` remains in the build for SGSIM and block
kriging support, but no longer supplies the PP/interbed fallback for model/data
SK. Its block-kriging implementation is retained, although libkriging OK/SK
currently reject block mode explicitly.

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
