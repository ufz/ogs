# OGS insitu-visualization with ParaView Catalyst

## Getting started

Requirements:

- ParaView insitu-build, at least commit [056de649](https://gitlab.kitware.com/paraview/paraview/commit/056de649320f52c8a14668ffa383d7361313a133), and therefore Conan is not supported

### Build ParaView

```bash
git clone --recursive https://gitlab.kitware.com/paraview/paraview.git
mkdir paraview_build paraview_install
cd paraview_build
cmake -G Ninja -DCMAKE_BUILD_TYPE=Release -DPARAVIEW_USE_PYTHON=ON \
  -DPARAVIEW_BUILD_EDITION=CATALYST -DCMAKE_INSTALL_PREFIX=../paraview_install \
  ../paraview
ninja install
```

### Build OGS

```bash
cmake -DCMAKE_BUILD_TYPE=Release -DOGS_INSITU=ON -DOGS_USE_CONAN=OFF \
  -DParaView_DIR=~/path/to/paraview_install/lib/cmake/paraview-5.8 \
  -DOGS_BUILD_PROCESSES=GroundwaterFlow ../ogs
```

**OR:** Build with [ogs-container-maker](https://github.com/ufz/ogs-container-maker):

```bash
python ogscm/cli.py --compiler_version 9 --ogs bilke/ogs@insitu-refactor --cmake_args ' -DOGS_BUILD_PROCESSES=GroundwaterFlow' --pm system --insitu -B -C -R
```

### Run benchmark

```
bin/ogs -o _out ../ogs/Tests/Data/Elliptic/cube_1x1x1_GroundWaterFlow/cube_1e1.prj
```

Open generated `cube_1e1_*.pvtp` in ParaView.

## How it works

See section `<insitu>` in `Elliptic/cube_1x1x1_GroundWaterFlow/cube_1e1.prj`. It defines Python scripts with visualization pipelines which are executed after every time step.

These python scripts can be generated with ParaView:

- Load a typical dataset representative for your expected simulation output
- Rename dataset in the *Pipeline Browser* to `input`
- Setup the filter pipeline
- In the menu click *Catalyst / Define Exports* which opens the *Export Inspector*
- In the *Export Inspector* under *Data Extracts*:
  - Select the filter you want to write out
  - Choose the data format, e.g. *XMLPPolyDataWriter*
  - Click the check-box next to the data format drop-down
- In the menu click *Catalyst / Export Catalyst Script*

----

TODO:

- Test parallel benchmark
- Live connection not working: https://gitlab.kitware.com/paraview/paraview/issues/19613
- Check https://gitlab.kitware.com/paraview/paraview/tree/master/Examples/Catalyst/CxxMappedDataArrayExample for more syntax changes
