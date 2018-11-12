#!/usr/bin/python

import vtk
from vtk.util.numpy_support import numpy_to_vtk
import numpy as np

reader = vtk.vtkXMLUnstructuredGridReader()
reader.SetFileName("cube_1x1x1_hex_1e0.vtu")
reader.Update()
grid = reader.GetOutput()

nnodes = grid.GetNumberOfPoints()

ts = np.arange(0, 11)
T0 = 273.15
Ts = np.linspace(273.15, 313.15, 9)
Ts = np.hstack(([ T0 ], Ts, [ Ts[-1] ]))
assert len(ts) == len(Ts)

E = 5e9
nu = 0.2
alpha = 1e-5
I = np.matrix([1, 1, 1, 0, 0, 0]).T

D = 0.0  # just a dummy value because that part of C is not needed!

C = np.matrix([
    [1-nu,  nu,    nu,    0,  0,  0],
    [nu,    1-nu,  nu,    0,  0,  0],
    [nu,    nu,    1-nu,  0,  0,  0],
    [0,     0,     0,     D,  0,  0],
    [0,     0,     0,     0,  D,  0],
    [0,     0,     0,     0,  0,  D]]) * (E / ( 1+nu ) / (1-2*nu))

def eps_T(T):
    return I * (alpha * (T - T0))

def sigma(T):
    return - C * eps_T(T)

for t, T in zip(ts, Ts):
    sigma_ = np.tile(sigma(T).T, (nnodes, 1))

    sigma_vtk = numpy_to_vtk(sigma_, 1)
    sigma_vtk.SetName(f"sigma_{t:g}")

    grid.GetPointData().AddArray(sigma_vtk)

writer = vtk.vtkXMLUnstructuredGridWriter()
writer.SetFileName("cthex_ref.vtu")
writer.SetInputData(grid)
writer.Write()
