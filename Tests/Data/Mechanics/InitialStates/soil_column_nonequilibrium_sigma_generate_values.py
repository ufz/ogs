#!/usr/bin/env python
from vtk import *


def createStressArray(points):
    sigma = vtkDoubleArray()
    sigma.SetNumberOfComponents(4)
    sigma.SetNumberOfTuples(points.GetNumberOfPoints())
    for i in range(points.GetNumberOfPoints()):
        y = points.GetPoint(i)[1]
        sigma.SetTuple4(
            i, 2200 * 9.81 * 0.8 * (y - 100), -0.5 * 1e6 + 2200 * 9.81 * (y - 100), 0, 0
        )
        # sigma.SetTuple4(i, 4*i, 4*i + 1, 4*i + 2, 4*i + 3) # for debugging
    sigma.SetName("nonequilibrium_stress")
    return sigma


def createStressArrayForCells(mesh):
    cell_centers = vtkCellCenters()
    cell_centers.SetInputData(mesh)
    cell_centers.Update()
    return createStressArray(cell_centers.GetOutput())


def createStressArrayForNodes(mesh):
    return createStressArray(mesh.GetPoints())


def writeDataToFile(mesh, filename):
    w = vtkXMLUnstructuredGridWriter()
    w.SetFileName(filename)
    w.SetInputData(mesh)
    w.Write()


def addStressArraysToMesh(mesh):
    mesh.GetPointData().AddArray(createStressArrayForNodes(mesh))
    mesh.GetCellData().AddArray(createStressArrayForCells(mesh))


r = vtkXMLUnstructuredGridReader()
r.SetFileName("soil_column.vtu")
r.Update()
m = r.GetOutput()
addStressArraysToMesh(m)
writeDataToFile(m, "soil_column_nonequilibrium_sigma.vtu")
