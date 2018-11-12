#!/usr/bin/python

import vtk
from vtk.util.numpy_support import numpy_to_vtk
import numpy as np

reader = vtk.vtkXMLUnstructuredGridReader()
reader.SetFileName("cube_1x1x1_hex_1e0.vtu")
reader.Update()
grid = reader.GetOutput()

nnodes = grid.GetNumberOfPoints()

ref_A = np.loadtxt("bdt.res")
ts = ref_A[:,0]
epss = ref_A[:,1:(1+6)]
sigmas = ref_A[:,(1+6):(1+12)]

for t, eps, sigma in zip(ts, epss, sigmas):
    # assert that we don't have to bother with permuting the shear
    # components from MFront ordering to VTK ordering
    assert all(eps[3:] == 0)
    assert all(sigma[3:] == 0)

    eps_ = np.tile(eps, (nnodes, 1))
    sigma_ = np.tile(sigma, (nnodes, 1))

    # print(eps_)
    # print(sigma_)

    eps_vtk = numpy_to_vtk(eps_, 1)
    eps_vtk.SetName(f"epsilon_{t:g}")

    sigma_vtk = numpy_to_vtk(sigma_, 1)
    sigma_vtk.SetName(f"sigma_{t:g}")

    grid.GetPointData().AddArray(eps_vtk)
    grid.GetPointData().AddArray(sigma_vtk)

writer = vtk.vtkXMLUnstructuredGridWriter()
writer.SetFileName("bdt_ref.vtu")
writer.SetInputData(grid)
writer.Write()
