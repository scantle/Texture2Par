# Res2Par integration with libkriging

## Objective

Port the contained `dev_res2par` preprocessing workflow onto the authoritative
`master` implementation without changing ordinary Texture2Par behavior when a
`PREPROC` block is absent.

The supported workflow interpolates log-resistivity to the model grid, applies
an optional per-cell multiple of the returned kriging standard deviation, and
converts the adjusted resistivity to normalized texture-class likelihoods.

## Design boundaries

- Preserve libkriging's interpolation mathematics. Its returned conditional
  variance is consumed by the preprocessor after interpolation; any library
  edits are limited to replacing large temporary/automatic arrays with heap or
  scalar-loop equivalents to prevent Windows stack exhaustion.
- Activate all new behavior only when a `PREPROC` block is present.
- Use a dedicated single-primary-variable preprocessing path rather than
  changing normal texture/covariate bookkeeping.
- Keep model-sized work arrays allocatable so they are heap allocated; avoid
  automatic arrays whose size depends on cells, observations, or classes.
- Preserve the current `master` behavior and version when preprocessing is off.
- Keep legacy one-depth-column conversion in test/model preparation rather than
  adding a second production dataset parser.

## Implementation

1. Add `src/PreProc.f90` with:
   - the base preprocessor type and `RES2TEXTURE` implementation;
   - validated scale/shape and multiplier file readers;
   - standard-deviation multiplier application using
     `sqrt(max(variance, 0.0))`;
   - NODATA/inactive-cell handling;
   - normalized lognormal likelihood conversion.
2. Add an internal `PREPROC_ON` option and own the preprocessor from `t_T2P`.
3. Parse the `PREPROC` input block and require it after `CLASSES` but before
   `VARIOGRAMS` and dataset interpretation.
4. Resolve preprocessing dataset and variogram class names through the
   preprocessor, including automatic and explicit dataset column mappings.
5. Adapt model interpolation so a preprocessing run interpolates one raw
   variable, retains libkriging variance, applies the preprocessor, then
   continues through the existing texture-to-parameter workflow.
6. Add the module to the makefile and Visual Studio builds in dependency order.
7. Remove observation-sized vector temporaries and the rotation routine's
   point-sized automatic array where the SVIHM run exposes stack exhaustion.

## Verification

- Unit/focused executable tests:
  - parsing and activation of `PREPROC`;
  - no-preprocessor behavior remains unchanged;
  - zero and nonzero standard-deviation multipliers;
  - normalized texture likelihoods;
  - malformed/incomplete scale and multiplier inputs;
  - NODATA and inactive cells.
- Run the existing pytest suite with the newly built executable.
- Prepare an isolated modern-format copy of the SVIHM legacy dataset:
  - first interval top is zero;
  - subsequent tops are the previous interval bottom;
  - legacy `Depth` becomes `BotDepth`;
  - omit the two zero-thickness records that the legacy reader effectively
    skipped.
- Compare SVIHM in stages:
  1. raw interpolated `LOGRHO` and kriging variance;
  2. multiplier-adjusted log-resistivity;
  3. five texture likelihood fields;
  4. final `KH`, `KV`, `SS`, and `SY` fields.

Any material raw kriging discrepancy will be investigated before adjusting
comparison tolerances, particularly for 3D block support, anisotropic search,
and measurement-error variance.
