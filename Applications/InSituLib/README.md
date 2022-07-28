# OGS insitu-visualization with ParaView Catalyst

## Getting started

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
cmake -DCMAKE_BUILD_TYPE=Release -DOGS_USE_INSITU=ON \
  -DParaView_DIR=~/path/to/paraview_install/lib/cmake/paraview-5.10 \
  -DOGS_BUILD_PROCESSES=SteadyStateDiffusion ../ogs
```

**OR:** Build with [ogs-container-maker](https://github.com/ufz/ogs-container-maker):

```bash
ogscm compiler.py ogs.py --cmake_args ' -DOGS_BUILD_PROCESSES=SteadyStateDiffusion' --insitu -B -C -R
```

### Run benchmark

```bash
bin/ogs -o _out ../ogs/Tests/Data/EllipticPETSc/square_1e1_neumann.prj
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
