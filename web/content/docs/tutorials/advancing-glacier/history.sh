#!/usr/bin/env bash

python3 -m venv .venv
source .venv/bin/activate
pip install --pre --index-url https://gitlab.opengeosys.org/api/v4/projects/120/packages/pypi/simple ogs gmsh ogstools

#create the mesh via python control for gmsh
python mesh_basin.py

#transform the mesh and extract boundary meshes
msh2vtu mesh_basin.msh --ogs --rdcd

#run OpenGeoSys (with the debug level information)
ogs -l debug OGSinput_basin.prj

#do the postprocessing with ParaView, disabled for CI
# paraview OGSoutput_basin0.pvd

deactivate
