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

### Recording executable timings

Set `T2P_TIMINGS_FILE` to record the wall time of each Texture2Par invocation,
excluding pytest setup, Python reference calculations, and assertions:

```powershell
$env:T2P_EXE = (Resolve-Path bin\Texture2Par.exe).Path
$env:T2P_TIMINGS_FILE = (Join-Path (Get-Location) timings.json)
pytest -q -p no:cacheprovider
```

The JSON output identifies the pytest node, invocation number, elapsed seconds,
and executable return code. Tests that invoke Texture2Par more than once produce
one record per invocation.

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
| `TestBatchModelIDW` | batch-backed model/data IDW exact hits, missing-value masks, and saved grid weights |
| `TestLibKrigingOrdinary` | default OK data, pilot-point, and interbed routing; OK_LIB alias; rejection of removed OK_LEGACY; explicit libkriging/GSLIB anisotropy and anisotropic-search coverage; NNEAR; nested structures; supported variogram mappings; duplicate and unsupported-feature diagnostics |
| `TestBlockKriging` | structured-MODFLOW block OK estimate/variance against an independent weighted-support reference; shared OK_LIB/SK routing; diagnostics for other models, IDW, and SGSIM |
| `TestLibKrigingSimple` | default SK routing through the shared libkriging backend; global-mean 2D/3D references; simple-cokriging secondary influence and per-variable means; SK/OK distinction; independent libkriging-OK pilot-point route; unsupported-feature diagnostics |
| `TestCorrectWeights` | `CORRECT_WEIGHTS` clips negative univariate kriging weights, renormalizes them to one, and changes estimates to an independent corrected-weight reference |
| `TestWriteMatrix` | `WRITE_MATRIX` emits libkriging neighbour, kriging-matrix, and right-hand-side CSV diagnostics |
| `TestCoincidentLocations` | duplicate XY/XYZ conditioning points are fatal; exact target hits remain valid |
| `TestSGS` | libkriging-backed OK/SK/CK 2-D/3-D point SGSIM; finite output; seed reproducibility; exact/soft conditioning variance; category masking; simple co-simulation means and secondary influence; unsupported-mode diagnostics; no legacy runtime allocation |
| `TestWeightCaching` | libkriging grid/PP saved-weight write/read, changed values, variance, CK, SK, block targets, and missing files |
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

### Libkriging feature matrix

| Feature | Current support |
|---|---|
| OK | Production libkriging path for model/data point support in 2-D/3-D and deterministic PP/interbed interpolation; `OK_LIB` is an alias. |
| SK | Deterministic model/data point support in 2-D/3-D and structured blocks, including simple cokriging with a known mean per variable; PP/interbed SK is unsupported because those workflows define no SK mean. |
| CK | Deterministic model/data point support in 2-D/3-D and structured blocks, including multiple secondary variables. |
| CK mean modes | `LIBKRIGING` and `LEGACY_STANDARDIZED`; both use the libkriging CK solve. |
| SGSIM OK | One point-support realization in 2-D/3-D, including categories/zones and observation variance. |
| SGSIM SK | One point-support realization in 2-D/3-D, including categories/zones, observation variance, and simple co-simulation with one known global mean per variable. |
| SGSIM CK | One primary-variable point-support realization in 2-D/3-D with CK mean modes, categories/zones, and observation variance. |
| `BLOCK_MODE` | Deterministic OK/`OK_LIB`/SK/CK for structured MODFLOW cells through libkriging Gaussian block support. |
| Kriging variance | End-to-end for deterministic OK/SK/CK point and block solves; grid factor files preserve it. |
| Observation variance | End-to-end for deterministic OK/SK/CK and SGSIM; deterministic factor files bind the write-run values. |
| Corrected weights | `CORRECT_WEIGHTS` uses libkriging's univariate negative-weight clipping and normalization for model/data, PP, and interbed solves; libkriging deliberately leaves multivariable CK weights unchanged. |
| Matrix diagnostics | `WRITE_MATRIX` writes libkriging `data_*.csv`, `matA_*.csv`, and `rhsB_*.csv` diagnostic files. |
| Saved grid weights | Deterministic OK/`OK_LIB`/SK/CK, including structured block targets, through libkriging factor files. |
| Saved PP weights | Deterministic PP OK through libkriging factor files. |
| Saved IB weights | Deterministic interbed OK through libkriging factor files. |
| IDW paths | Model/data, PP, and interbed use the shared batch interface. IDW does not support blocks or SGSIM. |
| Unsupported/future | Multiple SGSIM realizations; block SGSIM; PP/interbed SGSIM; saved SGSIM weights/paths; custom/unstructured blocks; sector search; asymmetric cross-variograms; factor-file versioning and geometry fingerprints. |

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
compatibility alias for that same path. `OK_LEGACY` has been removed. Any
anisotropic default-OK test must state which convention it
validates and must not reuse legacy
anisotropic expectations without explicitly documenting the expected
difference. The two sets of anisotropic tests must not be merged into one
parity expectation. The default-OK cutover intentionally changes anisotropic
results where the historical rotation convention differs.

`ANISOTROPIC_SEARCH` controls only the libkriging OK/SK neighbor-selection metric;
covariance evaluation continues to use every structure's own variogram
anisotropy regardless of that flag. Search currently uses the first variogram
structure's major/minor range ratios and angles. Sector search remains future
work and is not exposed by Texture2Par.

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

Saved deterministic kriging weights use libkriging's native factor-file system.
The existing `WRITE_GRID_WEIGHTS`/`READ_GRID_WEIGHTS`,
`WRITE_PP_WEIGHTS`/`READ_PP_WEIGHTS`, and
`WRITE_IB_WEIGHTS`/`READ_IB_WEIGHTS` options remain the user interface. Model
OK, `OK_LIB`, SK, CK, and structured `BLOCK_MODE` are supported, as are PP OK
and interbed OK. `READ_GRID_WEIGHTS` can be combined with `WRITE_VARIANCE`:
libkriging stores and restores the kriging variance along with neighbour
indices and weights. Observation variance is therefore baked into a factor
file and must not change when that file is reused. Observation values may
change; the loaded weights are applied to the current values.

Files retain the historical project/type stems but add a category suffix
because each category is a separate libkriging solve:
`<project>_grdweights_arrNN_obsN_catN.wts`,
`<project>_ppweights_arrNN_catN.wts`, and
`<project>_ibweights_layNN_catN.wts`. PP and interbed value columns share the
single category file and reuse its loaded in-memory weights. Native files are
formatted text whether or not the historical optional `TEXT` suffix is used.
Pre-libkriging `.wts` contents are not compatible and should be regenerated.
Target and observation counts and ordering, categories, search limits, support
mode, variograms, observation variance, and CK/SK configuration must match the
write run; the native header diagnoses dimension mismatches but cannot verify
all of that semantic geometry.

`SGSIM` is a simulation mode of the shared `t_kriging_interpolator` backend, not a separate
interpolator. `OK + SGSIM`, the equivalent `OK_LIB + SGSIM`, `SK + SGSIM`, and CK + SGSIM
reuse the same observation/target/category packing, variogram mapping, search
configuration, observation variance, solve, and scatter path as their
deterministic kriging modes. SK passes `unbias=0` and the existing known global
mean for every packed primary/secondary variable. Simple co-simulation uses
mean-centered residuals for both hard data and previously simulated nodes. CK
passes its existing multivariable observations,
auto/cross variograms, searches, and mean-mode configuration to libkriging's
native joint simulation and returns only the primary variable. The adapter
initializes libkriging with one simulation and the
Texture2Par `SEED`, supplies an in-memory randomized path and Gaussian samples,
calls `set_sim`, and returns the first simulated-value plane through the normal
result structure. Both `INTERP_DIM 2` and `3` and category/HSU solves are
supported. The current input flag and output arrays represent exactly one
realization.

SGSIM currently applies only to model/data point-support OK, SK, and CK
interpolation. IDW, `BLOCK_MODE`, saved weights or simulation paths, and
separate PP/interbed simulation are unsupported and fail clearly where those
combinations can be requested; no saved-path input option is exposed and the
native path remains in memory. Deterministic PP/interbed interpolation remains
independently controlled by `PP_INTERP_METHOD` and `IB_INTERP_METHOD` during an
SGSIM run. Observation variance is passed to libkriging, so zero variance gives
exact conditioning at coincident targets while positive variance provides soft
conditioning. Realizations are not expected to match the retired Texture2Par
solver bit-for-bit.

Default `INTERP_METHOD SK` now uses the same `t_kriging_interpolator` backend as OK. The
backend mode selects `unbias=1` for OK and `unbias=0` for SK; SK passes the
dataset/class global mean to libkriging's `set_obs`. That mean is computed by
the existing dataset readers across all layers/categories (and is
thickness-weighted for 3D interval data), not recomputed per category or
target. SK uses the libkriging search path and `NNEAR`; the removed
`t_globkrige` instead assembled all observations in the current category.

SK follows libkriging's GSLIB-style anisotropy convention, not Texture2Par's
historical rotation behavior. Because PP/interbed methods default independently
to OK, `INTERP_METHOD SK` routes both workflows through libkriging ordinary
kriging. `INTERP_ORDER AFTER` remains rejected for SK rather than silently
falling back because the generated parameter data never had an initialized
legacy SK mean.

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

Libkriging CK follows libkriging/GSLIB-style covariance and search anisotropy.
Simple cokriging uses each primary/secondary dataset's global mean and the same
multivariable packing, covariance, and search path as ordinary CK. Asymmetric
cross-variograms and sector search remain unsupported. `CK_MEAN_MODE` applies
only to ordinary CK; simple cokriging always uses its per-variable known means.
Saved CK weights are supported for deterministic runs; retain the same
`CK_MEAN_MODE` when reusing them. SGSIM support and limits are described above.

`BLOCK_MODE` historically used 4-point Gaussian quadrature per spatial axis
(16 support points in 2-D and 64 in 3-D) for structured MODFLOW cells. The
legacy OK/SK/CK routines averaged observation-to-support covariance for the
estimate and support-to-support covariance for the variance. It now maps to
libkriging's native Gaussian block API (`block_type=-4`): Texture2Par supplies
only block centres and half-dimensions, while libkriging generates the support
points and weights and performs both covariance averages. Deterministic `OK`,
`OK_LIB`, `SK`, and CK share the existing libkriging solve path in both 2-D and
3-D.

Texture2Par does not currently expose custom block support points or weights.
That may be added later for irregular grids such as IWFM. For now, block mode
remains restricted to structured MODFLOW because other flow-model readers do
not define cell dimensions. IDW and SGSIM are rejected explicitly rather than
performing point interpolation under a block label. Native libkriging currently
integrates hard-data-to-block covariance and each target block's conditional
variance, but previously simulated blocks re-enter the sequential system at
their centre coordinates rather than through block-to-block covariance. Block
SGSIM therefore remains disabled until that joint covariance is consistently
integrated. Pilot-point and interbed
interpolation remain point-support workflows even when model/data block mode is
active.

The old block covariance calculations and support-point storage have been
removed. The runtime no longer allocates `t_sgsim` or uses `t_treekrige` for
SGSIM. The legacy `sgsim.f90`, `kriging.f90`, and `kriging_vecback.f90` sources
have been retired and are absent from the makefile production builds;
production interpolation is libkriging-backed. The `kriging.F90` compiled under
`src/Libraries/libkriging` is libkriging's active implementation, not the
retired top-level source. Stale Visual Studio source entries remain pending a
separate project-file update.

Saved weights remain unsupported for SGSIM and saved simulation paths. Custom
or unstructured block support is not exposed. Model/data, PP, and interbed IDW
use the shared batch request and the dense factor-file format in `weights.f90`;
deterministic kriging factor files do not use that format.

Future factor-file hardening should add an explicit format version and geometry
fingerprint. This is separate from the currently enforced dimensions and
configuration metadata.

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
