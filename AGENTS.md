# Agent instructions for OpenGeoSys-6

**Project**: Scientific THMC simulation framework in modern C++. See [README.md](./README.md) for overview.

## Mandatory code standards

**C++**: C++23 standard (required). Use `ranges-v3` (not `std::ranges`), Eigen (not raw loops). Follow [style guide](https://ufz.github.io/styleguide/cppguide.xml).

**CMake**: 3.31+, target-based commands, no global variables.

**Python**: black-based formatting. ruff check linter.

**Documentation**: British English spelling (colour, behaviour).

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

- Build configuration via cmake presets, e.g. `release`, `release-petsc` for parallelised runs, build directory then in `../build/[preset]`
- Build with `ninja`
- Unit tests: `Tests/{LibName}/` (Google Test), `test`-target
- Integration tests: `Tests/Data/{ProcessName}/` (`.prj` files with reference outputs or `.xml`-patch files referencing `.prj`-files), typical test command in build directory: `ctest -LE large --output-on-failure`, full logs in [build-dir]/logs
- Always run ctests from release build.
- Python-based (`.py`) notebook tests (require `release-all` preset) also based on ctest, prefixed with `nb-`
- Check `.clang-format`, `.clang-tidy` for linting rules
- Use `pre-commit` to lint changes

## Documentation

- Doxygen C++ documentation, `doc`-target
- Hugo-based documentation web site, `web`-target, `preview-web`-target for local preview
