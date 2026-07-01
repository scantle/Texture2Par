# Texture2Par agent instructions

## Project context

Texture2Par is a modern Fortran codebase. The current migration goal is to replace the legacy Texture2Par kriging/SGSIM implementation with native Fortran libkriging while preserving Texture2Par input/output behavior through staged, tested changes.

Use only the local repository state. Do not fetch Texture2Par source files from GitHub unless explicitly asked.

## Editing rules

- Do not push to GitHub unless explicitly asked.
- Do not modify Visual Studio project files unless explicitly asked.
- Preserve existing line endings in every file touched. Some existing files use CRLF. Do not introduce mixed CRLF/LF endings or noisy whole-file rewrites.
- Before handoff, inspect the final diff. If line-ending churn was introduced, normalize only the files you touched back to their original line-ending style.
- Prefer small, compiling, tested changes.
- If a requested change becomes much larger than expected or requires speculative architecture decisions, stop and ask.

## Build rules

The Intel Fortran makefile build is under `build_ifx/makefile`.

Typical release rebuild command from `build_ifx`:

```powershell
cmd /c 'call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64 >nul && make -f makefile CFG=release rebuild'
```

Do not modify Visual Studio project files as part of makefile-oriented tasks.

## Testing rules

The pytest environment is `t2p_tests`.

- Previous known-good executable: `.\bin\Texture2Par_prev.exe`
- Candidate executable: `.\bin\Texture2Par.exe`
- `T2P_EXE` is authoritative. Always confirm pytest reports the intended executable path.

Typical candidate full-suite command from the repository root:

```powershell
$env:T2P_EXE=(Resolve-Path '.\bin\Texture2Par.exe').Path; conda run -n t2p_tests pytest --basetemp=.pytest_tmp -x --tb=long
```

Typical previous-executable baseline command:

```powershell
$env:T2P_EXE=(Resolve-Path '.\bin\Texture2Par_prev.exe').Path; conda run -n t2p_tests pytest --basetemp=.pytest_tmp -x --tb=long
```

Because the legacy KD-tree has temporarily been replaced by brute-force search, the full pytest suite currently takes about 9 minutes 20 seconds. Use focused tests during development and a timeout of at least 20 minutes for full-suite runs.

A `.pytest_tmp` cleanup/cache warning may occur due to Windows ACL issues. Do not treat it as a failure if tests otherwise pass.

## Migration guardrails

- Keep existing `OK`, `SK`, `IDW`, and SGSIM behavior unchanged unless the task explicitly says to change them.
- `OK_LIB` is the experimental libkriging ordinary-kriging path.
- Do not reproduce the old OK/SK duplication in new libkriging code. libkriging uses one `t_kriging` type with OK/SK as configuration.
- Do not silently switch OK to SK for efficiency. Preserve the requested kriging mode.
- Texture2Par's legacy anisotropy/dip convention differs from libkriging/GSLIB/scipy-style rotation. Do not force libkriging to match the legacy convention unless explicitly asked.
- Pilot-point values are stored as `values(point, parameter, class_id)`, not like layerdata `values(class, point)`.
- In `INTERP_DIM 3`, observations are treated as 3D support points and target Z values correspond to model-layer center elevations. Do not accidentally layer-average 3D data like 2D layer interpolation.
- Duplicate/co-located conditioning observations are fatal. A target coincident with a single observation is allowed.
