# Migrating cokriging models from Texture2Par 2.0.x to 2.1

Texture2Par 2.1 uses the KrigeKit/libkriging engine for kriging and cokriging. This is a methodological and numerical update, not a bit-for-bit replacement of the 2.0.x interpolation implementation. Users should expect cokriging results to change, even when input datasets and variograms are unchanged.

## Summary

For new models, the 2.1 cokriging formulation is recommended. It preserves the multivariable constraints of the cokriging system and does not apply the legacy negative-weight clipping and joint renormalization used by 2.0.x when `CORRECT_WEIGHTS` was enabled.

For an existing calibrated model:

1. Archive the 2.0.x executable, inputs, outputs, and log needed to reproduce the historical model.
2. Run 2.1 using its default `CK_MEAN_MODE LIBKRIGING` behavior.
3. Compare texture and derived hydraulic-property fields spatially, not only with domain-wide summary statistics.
4. Run `CK_MEAN_MODE LEGACY_STANDARDIZED` as a sensitivity test when the primary and secondary datasets have different global means.
5. Revisit calibration if the 2.0.x texture or parameter fields were part of a calibrated groundwater model.

`LEGACY_STANDARDIZED` is a mean-treatment compatibility setting; it is not a complete 2.0.x numerical compatibility mode.

## Cokriging weight constraints

In ordinary cokriging with one primary and one secondary variable, the standard unbiasedness constraints require the primary-data weights to sum to one and the secondary-data weights to sum to zero. Individual primary or secondary weights may legitimately be negative.

When `CORRECT_WEIGHTS` was enabled, the 2.0.x implementation could:

1. replace negative primary and secondary observation weights with zero; and
2. renormalize all remaining observation weights together to sum to one.

That operation can be useful as an empirical smoothing rule, but it no longer represents the constrained cokriging solution: it changes the separate weight-sum constraints and can make the prediction depend on an arbitrary offset in the secondary variable.

Texture2Par 2.1 therefore does not apply `CORRECT_WEIGHTS` to multivariable cokriging. The program writes a warning when `CORRECT_WEIGHTS` is requested for a cokriging model. The option remains available for univariate kriging, where clipping and renormalizing one set of weights has a different interpretation.

This change is commonly the largest reason that a 2.1 cokriging result differs from its 2.0.x counterpart. It is intentional and is considered the more rigorous default formulation.

## Cokriging mean modes

Texture2Par 2.1 recognizes the following main-input options:

```text
CK_MEAN_MODE LIBKRIGING
CK_MEAN_MODE LEGACY_STANDARDIZED
```

### `LIBKRIGING`

This is the default and recommended mode for new work. It uses the KrigeKit/libkriging mean treatment together with the multivariable cokriging constraints.

### `LEGACY_STANDARDIZED`

This mode shifts each secondary dataset so that its global mean is aligned with the primary dataset mean before interpolation, following the historical Texture2Par convention. It can be useful when:

- comparing a migrated model with archived 2.0.x results;
- testing sensitivity to a global offset between borehole and geophysical datasets; or
- representing a documented measurement-scale or support-scale correction.

It should not be selected solely because it makes a new result resemble an old result. The global mean difference may represent real sampling, spatial-coverage, or support differences rather than a bias that should be removed.

Because `LEGACY_STANDARDIZED` still uses the 2.1 constrained cokriging weights, it will not exactly reproduce a 2.0.x run that used legacy corrected weights.

## Data locations outside the model mesh

Version 2.1 uses a newer deterministic grid-search implementation. Observations inside the mesh retain a geometrically defined element or cell. An observation outside the mesh but permitted by `MAX_OUTSIDE_DIST` can be equidistant from multiple boundary elements, making its assignment geometrically ambiguous.

Models that use outside-mesh observations—particularly vertically averaged AEM data, for which the selected element supplies model-layer elevations—should review their locations and resulting layer averages during migration. Prefer observations inside the model domain or use a documented preprocessing rule rather than relying on search order to resolve an exact boundary tie.

## Missing well identifiers

Version 2.1 initializes output location identifiers for missing or discontinuous well numbers. This prevents an uninitialized element number from appearing in layer-averaged diagnostic output for a missing well. The correction affects diagnostic bookkeeping rather than the interpolation result for valid observations.

## Recommended interpretation of changed results

A difference from 2.0.x is not, by itself, evidence that the 2.1 result is incorrect. Evaluate the migrated model using:

- cross-validation or held-out observations;
- bias and error statistics appropriate to the modeled variables;
- maps and vertical sections of residuals and parameter changes;
- groundwater-model calibration and predictive behavior; and
- sensitivity to the two available cokriging mean modes.

For final reporting, record the Texture2Par version, `CK_MEAN_MODE`, interpolation method, variogram file, and whether `CORRECT_WEIGHTS` was requested.

## References

- [Texture2Par method paper](https://doi.org/10.1016/j.envsoft.2025.106372)
- [DWR Basin Characterization Tools](https://data.ca.gov/@california-department-of-water-resources/dwrs-basin-characterization-tools-bctools)
- [KrigeKit repository](https://github.com/ougx/krigekit)
- [KrigeKit v0.3.2 preprint](https://doi.org/10.5194/egusphere-2026-5708)
- [Probabilistic multi-texture AEM application](https://doi.org/10.1029/2025WR042695)
