# Texture2Par

Texture2Par is a Fortran groundwater-model preprocessor for constructing three-dimensional sediment-texture models and translating texture and aquifer-property data into parameter fields for MODFLOW-family models and the Integrated Water Flow Model (IWFM). It supports lithologic logs, airborne electromagnetic (AEM) interpretations, geologic or hydrostratigraphic zones, pilot points, kriging and cokriging, and empirical relationships between texture and hydraulic properties.

## Development and lineage

Texture2Par was developed by [S.S. Papadopulos & Associates, Inc. (SSP&A)](https://sspa.com/). It builds on an original tool and method developed by Tim Durbin and Claire Velayas of Timothy J. Durbin, Inc.

> [!IMPORTANT]
> This is an active development repository and contains changes made after the official California Department of Water Resources (DWR) release. For the official DWR distribution, manuals, examples, and training material, use [DWR's Basin Characterization Tools (BCTools)](https://data.ca.gov/@california-department-of-water-resources/dwrs-basin-characterization-tools-bctools).

## Capabilities

- Aggregate interval-based texture observations to model layers.
- Interpolate texture using ordinary, simple, block, and co-kriging workflows.
- Incorporate multiple texture classes and secondary datasets such as AEM-derived texture information.
- Convert texture fields to hydraulic conductivity and storage properties.
- Interpolate spatially varying aquifer parameters using pilot points.
- Calculate interbed thickness information used in subsidence parameterization.
- Write parameter inputs for IWFM and MODFLOW-family groundwater models.
- Run as a standalone preprocessor or within parameter-estimation and calibration workflows.

## Documentation and examples

The official BCTools distribution provides the current DWR user manuals, quick guide, example datasets, and training videos:

- [DWR BCTools downloads and documentation](https://data.ca.gov/@california-department-of-water-resources/dwrs-basin-characterization-tools-bctools)
- [Texture2Par method paper](https://doi.org/10.1016/j.envsoft.2025.106372)
- Repository examples: [`examples/`](examples/)
- Automated test-suite documentation: [`tests/README.md`](tests/README.md)

Because this repository includes post-release development, its executable behavior and available options may differ from the official DWR package. Use documentation from the same version as the executable being run.

## Building

The supported build paths currently target 64-bit Windows and Intel oneAPI Fortran:

### Visual Studio

Open `Texture2Par.sln` in Visual Studio 2022 with Intel oneAPI Fortran integration and build the `Release|x64` configuration.

### Command line

From an Intel oneAPI command prompt with GNU Make available:

```powershell
cd build_ifx
make -f makefile rebuild
```

If GNU Make is installed as `mingw32-make`, substitute that command for `make`. Use `CFG=debug` for a debug build. The makefile is currently verified only on Windows with Intel IFX.

Both build routes write `bin/Texture2Par.exe`. This repository does not currently provide a supported Linux or macOS build configuration.

### Intel OpenMP runtime DLL

The downloadable Windows release is built with OpenMP support and requires
`libiomp5md.dll`. This file is Intel's OpenMP runtime library; it provides the
thread creation, scheduling, and synchronization used when Texture2Par runs
deterministic libkriging solves in parallel. It is a program dependency, not a
Texture2Par input file or optional plugin.

Keep `libiomp5md.dll` in the same directory as `Texture2Par.exe`. Windows
searches the executable's directory for dependent DLLs, so this works even
when Texture2Par is launched from a different working directory. If the
executable is copied into a model or working directory, copy the DLL with it.
The DLL is still required when `NTHREADS` is 1 because Windows loads program
dependencies when the executable starts.

A non-OpenMP build, such as the Visual Studio `Release|x64` configuration,
does not require this DLL. The supported release package uses
`Release_OpenMP|x64`, however, so users should extract or copy both bundled
files. If Windows reports that `libiomp5md.dll` is missing, re-extract the
complete release ZIP or place its bundled copy beside the executable. Do not
download replacement DLLs from third-party DLL sites.

### Building a release

The release helper builds the 64-bit OpenMP configuration, copies the required
Intel OpenMP runtime, runs the full pytest suite against that exact executable,
and creates a versioned ZIP plus SHA-256 checksum under `build/release/`:

```powershell
.\release.bat 2.1.1
```

The version argument must match `version.f90`. Use `-SkipTests` only for local
packaging checks, never for a published release.

After installing and authenticating the [GitHub CLI](https://cli.github.com/),
the same command can tag the tested commit, push it, and create a GitHub Release:

```powershell
.\release.bat 2.1.1 -Publish -Draft
```

Omit `-Draft` when the release should be published immediately. Add
`-Prerelease` for a beta or release candidate. Publishing requires a clean
worktree; build and test outputs are ignored by Git.

## Testing

Install the Python test dependencies and run pytest against the compiled executable:

```powershell
py -m venv .venv-tests
.\.venv-tests\Scripts\python.exe -m pip install -r requirements-test.txt
$env:T2P_EXE = (Resolve-Path .\bin\Texture2Par.exe).Path
.\.venv-tests\Scripts\python.exe -m pytest
```

See [`tests/README.md`](tests/README.md) for test categories, binary selection, timing output, and targeted test commands.

## KrigeKit/libkriging

Texture2Par 2.1 uses a bundled and adapted copy of the Fortran kriging engine from [KrigeKit](https://github.com/ougx/krigekit). KrigeKit provides ordinary, simple, and co-kriging, geostatistical simulation, anisotropy, and related numerical capabilities. Its bundled license is available at [`src/Libraries/libkriging/LICENSE`](src/Libraries/libkriging/LICENSE).

For a description of KrigeKit's design and validation, see:

> Ou, G., Bedekar, V., Tonkin, M., Barth, G., and Karanovic, M. (2026). KrigeKit v0.3.2: a parallel Python–Fortran engine for variogram analysis, kriging, and geostatistical simulation. *EGUsphere* [preprint]. [https://doi.org/10.5194/egusphere-2026-5708](https://doi.org/10.5194/egusphere-2026-5708)

## Version 2.1 compatibility notice

Texture2Par 2.1 does **not** reproduce Texture2Par 2.0.x cokriging results exactly. Version 2.1 replaces the legacy cokriging implementation with the KrigeKit/libkriging engine and preserves the separate primary and secondary unbiasedness constraints of multivariable cokriging. Version 2.0.x could clip negative cokriging weights and jointly renormalize the remaining primary and secondary weights when `CORRECT_WEIGHTS` was enabled. That legacy correction changes the cokriging solution and is not applied to multivariable cokriging in 2.1.

Version 2.1 also provides two cokriging mean treatments:

- `CK_MEAN_MODE LIBKRIGING` is the default and is recommended for new work.
- `CK_MEAN_MODE LEGACY_STANDARDIZED` applies the historical global secondary-mean shift for sensitivity analysis and comparison with older models. It does not restore the legacy clipped-weight solution and therefore is not an exact 2.0.x compatibility mode.

Models calibrated using 2.0.x cokriging should be reviewed and, where appropriate, recalibrated before adopting 2.1 results. See the [Version 2.1 migration guide](documentation/version-2.1-migration.md) for the rationale, expected differences, and recommended workflow.

## Scientific references

The methods and original applications of Texture2Par are described in:

> Scantlebury, L., Bedekar, V., Tonkin, M. J., Karanovic, M., and Harter, T. (2025). Texture2Par: A texture-driven tool for estimating subsurface hydraulic properties. *Environmental Modelling & Software*, 187, 106372. [https://doi.org/10.1016/j.envsoft.2025.106372](https://doi.org/10.1016/j.envsoft.2025.106372)

A recent application integrating AEM and borehole information in a probabilistic multi-texture workflow is:

> Scantlebury, L., and Harter, T. (2026). From resistivity to hydraulic properties: Calibrating a groundwater flow model by integrating airborne electromagnetic and borehole data into a probabilistic multi-texture framework. *Water Resources Research*, 62, e2025WR042695. [https://doi.org/10.1029/2025WR042695](https://doi.org/10.1029/2025WR042695)

## Repository layout

- `src/` — Texture2Par and bundled-library source code
- `documentation/` — supplemental and version-specific documentation
- `examples/` — example inputs and model data
- `tests/` — end-to-end pytest regression suite
- `Texture2Par.sln` — Visual Studio solution
- `version.f90` — program version string
