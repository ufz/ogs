#!/usr/bin/vtkpython

import sys
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
import vtk

in_grid, out_grid, out_geom = sys.argv[1:]


def distribute_points_evenly(c2):
    assert np.sqrt(c2.shape[0]).is_integer()

    CELLS_PER_DIRECTION = int(np.sqrt(c2.shape[0])) - 1
    r2 = np.sqrt(c2[:,0]**2 + c2[:,1]**2)
    alpha2 = np.arctan2(c2[:,1], c2[:,0])

    bins = [ [] for _ in range(CELLS_PER_DIRECTION+1) ]
    nbin = np.round(r2 * CELLS_PER_DIRECTION).astype(int)
    for node, b in enumerate(nbin):
        bins[b].append(node)

    for b in bins:
        b.sort(key=lambda n: alpha2[n])
        # print(len(b))

    c3 = np.zeros_like(c2)

    for node, r in enumerate(r2):
        b = nbin[node]
        i = bins[b].index(node)
        if len(bins[b]) == 1:
            phi = 0.0
        else:
            phi = np.pi * 0.5 / (len(bins[b]) - 1) * i

        c3[node, 0] = r * np.cos(phi)
        c3[node, 1] = r * np.sin(phi)

    return c3


reader = vtk.vtkXMLUnstructuredGridReader()
reader.SetFileName(in_grid)
reader.Update()

grid = reader.GetOutput()
assert grid.GetBounds() == (0.0, 1.0, 0.0, 1.0, 0.0, 0.0)

coords = vtk_to_numpy(grid.GetPoints().GetData())

a = np.empty(coords.shape[0])

for i, (x, y, z) in enumerate(coords):
    if x > y:
        R = np.sqrt(1 + (y/x)**2)
        a[i] = 1.0 / R
    elif x < y:
        R = np.sqrt(1 + (x/y)**2)
        a[i] = 1.0 / R
    else:
        a[i] = np.sqrt(0.5)

# scale coordinates
new_coords = np.multiply(coords.T, a).T

if False:
    # If the input is a regular mesh, there is the possibility to make bins of
    # points with equal radius. For each radius those points can than be
    # distributed with even spacing at the quarter circles.
    new_coords = distribute_points_evenly(new_coords)

new_coords_vtk = numpy_to_vtk(new_coords, 1)
pts = vtk.vtkPoints()
pts.SetData(new_coords_vtk)
grid.SetPoints(pts)

writer = vtk.vtkXMLUnstructuredGridWriter()
writer.SetFileName(out_grid)
writer.SetInputData(grid)
writer.Write()

# extract boundary of unit circle
R_squared = new_coords[:,0]**2 + new_coords[:,1]**2
phi = np.arctan2(new_coords[:,1], new_coords[:,0])

idcs = np.where(abs(R_squared - 1) < 1e-8)[0]
idcs = idcs[np.argsort(phi[idcs])]  # sorted with ascending polar angle


with open(out_geom, "w") as fh:
    fh.write("""<?xml version="1.0" encoding="ISO-8859-1"?>
<OpenGeoSysGLI>
    <name>geom</name>
    <points>
        <point id="0" x="0" y="0" z="0" name="center"/>
        <point id="1" x="0" y="1" z="0" name="top"/>
        <point id="2" x="1" y="0" z="0"/>
""")

    for i, (x, y, z) in enumerate(new_coords[idcs]):
        fh.write('        <point id="{}" x="{}" y="{}" z="{}" />\n'.format(
            i+3, x, y, z))

    fh.write("""
    </points>

    <polylines>
        <polyline id="0" name="left">
            <pnt>0</pnt>
            <pnt>1</pnt>
        </polyline>
        <polyline id="1" name="bottom">
            <pnt>0</pnt>
            <pnt>2</pnt>
        </polyline>
        <polyline id="2" name="outer">\n""")

    for i in range(len(idcs)):
        fh.write("            <pnt>{}</pnt>\n".format(i+3))

    fh.write("""</polyline>
    </polylines>
</OpenGeoSysGLI>""")
