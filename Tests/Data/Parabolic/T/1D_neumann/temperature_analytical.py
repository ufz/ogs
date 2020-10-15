#!/usr/bin/env python
# Solution of heatequation in a semi-infinite domain.

from vtk import *
import numpy as np
from scipy.special import erfc

r = vtkXMLUnstructuredGridReader()
r.SetFileName("mesh.vtu")
r.Update()
m = r.GetOutput()
ps = m.GetPoints()
pd = m.GetPointData()

# PDE's coefficients:
lambda_coeff = 3.2

# alpha is the thermal diffusivity
alpha = lambda_coeff / (1000 * 2500)

T_inf = 273.15
q = 2


def temperature(x, t):
    return T_inf + 2 * q / lambda_coeff * (
        np.sqrt(alpha * t / np.pi) * np.exp(-(x ** 2) / (4 * alpha * t))
        - x / 2 * erfc(x / (2 * np.sqrt(alpha * t)))
    )


def addSolution(t):
    T = vtkDoubleArray()
    T.SetName("temperature_" + str(t) + "s")
    T.SetNumberOfComponents(1)
    T.SetNumberOfTuples(ps.GetNumberOfPoints())

    for i in range(ps.GetNumberOfPoints()):
        x = ps.GetPoint(i)[0]
        T.SetTuple1(i, temperature(x, t))

    pd.AddArray(T)


for t in 78125 * np.array([1, 3, 65, 405, 500]):
    addSolution(t)

w = vtkXMLUnstructuredGridWriter()
w.SetFileName("x.vtu")
w.SetInputData(m)
w.Update()
