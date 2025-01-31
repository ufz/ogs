#!/usr/bin/vtkpython

import xml.etree.ElementTree as ET
from pathlib import Path

import chamber as ch
import matplotlib.pyplot as plt
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

pvd_file = "out/piston_pcs_0.pvd"


### helpers ##############################################


def relpathfrom(origin, relpath):
    if relpath.is_absolute():
        return str(relpath)
    return str(origin / relpath)


def read_pvd_file(fn):
    try:
        path = fn.name
    except AttributeError:
        path = fn
    pathroot = Path(path).parent
    pvdtree = ET.parse(fn)
    node = pvdtree.getroot()
    if node.tag != "VTKFile":
        return None, None
    children = list(node)
    if len(children) != 1:
        return None, None
    node = children[0]
    if node.tag != "Collection":
        return None, None

    ts = []
    fs = []

    for child in node:
        if child.tag != "DataSet":
            return None, None
        ts.append(float(child.get("timestep")))
        fs.append(relpathfrom(pathroot, Path(child.get("file"))))

    return ts, fs


### helpers end ##########################################


ts, fns = read_pvd_file(pvd_file)

reader = vtk.vtkXMLUnstructuredGridReader()


loc_points = vtk.vtkPoints()
loc_points.InsertNextPoint([0.0, 0.0, 0.0])
loc = vtk.vtkPolyData()
loc.SetPoints(loc_points)

probe = vtk.vtkProbeFilter()
probe.SetSourceConnection(reader.GetOutputPort())
probe.SetInputData(loc)

uys = np.zeros(len(ts))
ps = np.zeros(len(ts))

for i, (t, fn) in enumerate(zip(ts, fns)):
    print("###### time", t)
    reader.SetFileName(fn)
    probe.Update()

    grid = probe.GetOutput()
    uy = vtk_to_numpy(grid.GetPointData().GetArray("displacement"))[0, 1]
    p = vtk_to_numpy(grid.GetPointData().GetArray("sigma"))[0, 1]

    uys[i] = uy
    ps[i] = -p

uys_ana = np.linspace(min(uys), max(uys), 100)
ps_ana = ch.p_chamber(uys_ana)


fig, ax = plt.subplots()
ax.scatter(uys[1:], ps[1:], label="ogs")
ax.plot(uys_ana, ps_ana, label="ref")
ax.legend()

ax.set_ylabel("pressure inside the chamber $p$ / Pa")
ax.set_xlabel("displacement of the piston $u_y$ / m")

fig.subplots_adjust(left=0.15, right=0.98)
fig.savefig("pressure-displacement.png")
plt.close(fig)
