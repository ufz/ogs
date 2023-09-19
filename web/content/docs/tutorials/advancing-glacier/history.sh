#create the mesh via python control for gmsh
python3 mesh_basin.py

#transform the mesh and extract boundary meshes
python3 msh2vtu.py mesh_basin.msh --ogs --rdcd

#run OpenGeoSys (with the debug level information)
myPATH2OGS=~/Forschung/gitprojects/OGS/build-release/bin
${myPATH2OGS}/ogs -l debug OGSinput_basin.prj

#do the postprocessing with ParaView
paraview OGSoutput_basin0.pvd
