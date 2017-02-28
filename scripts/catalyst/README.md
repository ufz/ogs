# Generate OGS Catalyst edition

Clone ParaView:

```bash
git clone https://gitlab.kitware.com/paraview/paraview.git
cd paraview
git submodule init
git submodule update
```

Go to ogs source directoy and run script `generate-ogs-edition.sh` with the path to the cloned ParaView directory as an argument:

```bash
cd ogs/scripts/catalyst
./generate-ogs-edition.sh ~/code/catalyst/paraview
```

This will generate the file `Catalyst-ogs-Base-Enable-Python-Essentials-Extras-ogs-Rendering-Base.tar.gz` in `[paraview-dir]/Catalyst/`

# Build the OGS Catalyst edition

- Untar the edition
- `cd` into it and create a `build`-directory
- Build with `cmake.sh`:

```bash
../cmake.sh .. [Optional CMake parameter]
make
```

# Use the OGS Catalyst edition

Enable `OGS_INSITU` and point to build directory where the edition was built:

```bash
cmake ../ogs -DOGS_INSITU=ON -DParaView_DIR=[path to edition build dir]
```
