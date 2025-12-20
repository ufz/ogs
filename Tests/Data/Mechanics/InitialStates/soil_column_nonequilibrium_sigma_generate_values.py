#!/usr/bin/env python
# SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
# SPDX-License-Identifier: BSD-3-Clause

import vtk


def createStressArray(points):
    sigma = vtk.vtkDoubleArray()
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
    cell_centers = vtk.vtkCellCenters()
    cell_centers.SetInputData(mesh)
    cell_centers.Update()
    return createStressArray(cell_centers.GetOutput())


def createStressArrayForNodes(mesh):
    return createStressArray(mesh.GetPoints())


def writeDataToFile(mesh, filename):
    w = vtk.vtkXMLUnstructuredGridWriter()
    w.SetFileName(filename)
    w.SetInputData(mesh)
    w.Write()


def addStressArraysToMesh(mesh):
    mesh.GetPointData().AddArray(createStressArrayForNodes(mesh))
    mesh.GetCellData().AddArray(createStressArrayForCells(mesh))


r = vtk.vtkXMLUnstructuredGridReader()
r.SetFileName("soil_column.vtu")
r.Update()
m = r.GetOutput()
addStressArraysToMesh(m)
writeDataToFile(m, "soil_column_nonequilibrium_sigma.vtu")
