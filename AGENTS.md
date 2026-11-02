# Agent instructions for OpenGeoSys-6

**Project**: Scientific THMC simulation framework in modern C++. See
[README.md](./README.md) for overview.

## Mandatory code standards

**C++**: C++23 standard (required). Use `ranges-v3` (not `std::ranges`), Eigen
(not raw loops). Follow
[style guide](https://ufz.github.io/styleguide/cppguide.xml).

**CMake**: 3.31+, target-based commands, no global variables.

**Python**: black-based formatting. ruff check linter.

**Documentation**: British English spelling (colour, behaviour).

**Markdown**: wrap prose at 80 columns; tables, long links, and code may exceed.

## Architecture layers

Foundation:  BaseLib → MathLib → NumLib
Geometry:    GeoLib, MeshLib, MeshGeoToolsLib, MeshToolsLib
Materials:   MaterialLib (MPL), ParameterLib
Processes:   ProcessLib (20+ process implementations)
Apps:        CLI, ApplicationsLib, FileIO, Utils, DataExplorer(Qt)

## Process implementation pattern (REQUIRED)

Every process must have:

1. `{Name}Process.h` - Inherits Process, manages assembly/timestepping
2. `{Name}ProcessData.h` - Material properties, parameters, solver configuration
3. `{Name}LocalAssembler.h` - Element-level assembly (M, K, b matrices)
4. `Create{Name}Process.h` - Factory from XML configuration

## Testing & validation

- Build configuration via cmake presets, e.g. `release`, `release-petsc` for
  parallelised runs, build directory then in `../build/[preset]`
- Build with `ninja`
- Unit tests: `Tests/{LibName}/` (Google Test), `test`-target
- Integration tests: `Tests/Data/{ProcessName}/` (`.prj` files with reference
  outputs or `.xml`-patch files referencing `.prj`-files), typical test command
  in build directory: `ctest -LE large --output-on-failure`, full logs in
  [build-dir]/logs
- Scope tests to the change: build the `testrunner` target and run
  `bin/testrunner --gtest_filter='<Suite>*'` for the touched libraries, plus
  `ctest -R <name>` for the benchmarks that exercise them. The full
  `ctest -LE large` takes hours; run it for changes to processes or shared
  infrastructure.
- Always run ctests from release build.
- Python-based (`.py`) notebook tests (require `release-all` preset) also based
  on ctest, prefixed with `nb-`
- Check `.clang-format`, `.clang-tidy` for linting rules
- Use `pre-commit` to lint changes

## Documentation

- Doxygen C++ documentation, `doc`-target
- Hugo-based documentation web site, `web`-target, `preview-web`-target for
  local preview

## Commit messages

Conventional commits: `type(scope): subject` — imperative, lower-case subject,
no trailing period; body explains *why*.

**Types**: `feat`, `fix`, `refactor`, `perf`, `test`, `docs`, `build`, `ci`,
`style` (formatting only), `chore`.

**Scope**: library or process abbreviation; sub-paths with `/`
(`NL/TimeStepping`, `MeL/IO`, `T/THM`); several scopes with `,`. The list below
is not exhaustive; extend it when needed.

- Libraries: `BL` BaseLib, `MaL` MathLib, `NL` NumLib, `GL` GeoLib, `MeL`
  MeshLib, `MGTL` MeshGeoToolsLib, `MTL` MeshToolsLib, `MatL` MaterialLib, `MPL`
  MaterialLib/MPL, `ParL` ParameterLib, `CL` ChemistryLib, `PL` ProcessLib,
  `AppL` ApplicationsLib, `App` Applications (`App/U` Utils, `App/DE`
  DataExplorer, `App/IO`), `T` Tests. Not `ML` (MeshLib or MathLib) nor `HC`
  (HeatConduction or ComponentTransport).
- Processes: `LF` LiquidFlow, `HT`, `CT` ComponentTransport, `RCT`
  RichardsComponentTransport, `RF` RichardsFlow, `RM` RichardsMechanics, `TRF`
  ThermoRichardsFlow, `TRM` ThermoRichardsMechanics, `SD` SmallDeformation, `LD`
  LargeDeformation, `HM`, `TM` ThermoMechanics, `THM`, `TH2M`, `LIE`, `BHE`
  HeatTransportBHE, `T2PF` TwoPhaseFlowWithPP, `TH2PP`
  ThermalTwoPhaseFlowWithPP, `PhF` PhaseField, `HMPhF` HMPhaseField, `SSD`
  SteadyStateDiffusion, `HeatConduction`, `WellboreSimulator`.
- ProcessLib parts: `PL/BC`, `PL/ST`, `PL/DS` DeactivatedSubdomain, `PL/Output`,
  `PL/TimeLoop`.
- Other: `CTest` test registration, `nb` notebooks, `py` Python, `web`, `doc`
  Documentation/ProjectFile (prj docs), `cmake`, `scr` scripts, `guix`,
  `pre-commit`, `agents` AGENTS.md and agent instructions, `wheel` Python
  wheels.

Examples: `fix(THM): ...`, `refactor(NL/DOF): ...`, `test(T/TRM): ...`.
