#!/usr/bin/vtkpython

import xml.etree.ElementTree as ET
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import vtk
from scipy.interpolate import interp1d
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy

pvd_file = "out/hertz_pcs_0.pvd"


R12 = 1.0
R = R12 / 2.0
nu_12 = 0.3
E_12 = 1.0

lambda_ = E_12 * nu_12 / (1 + nu_12) / (1 - 2 * nu_12)
mu = E_12 / 2.0 / (1 + nu_12)
G_12 = mu

kappa = 0.5 * G_12 / R / (1 - nu_12)
print("kappa:", kappa)

C = lambda_ * np.matrix([1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0]).reshape(
    4, 4
) + 2 * mu * np.identity(4)


def p_contact(r, r_contact):
    return kappa * np.sqrt(r_contact**2 - r**2)


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


def get_y_top(t):
    return 1.0 - 0.005 * t


ts, fns = read_pvd_file(pvd_file)

reader = vtk.vtkXMLUnstructuredGridReader()

destroyTopology = vtk.vtkShrinkFilter()
destroyTopology.SetShrinkFactor(1.0)
# destroyTopology.SetInputConnection(strainFilter.GetOutputPort())
# destroyTopology.SetInputConnection(reader.GetOutputPort())

# strainFilter = vtk.vtkCellDerivatives()
# strainFilter.SetVectorModeToPassVectors()
# strainFilter.SetTensorModeToComputeStrain()
# strainFilter.SetOutputPointsPrecision(vtk.vtkAlgorithm.DOUBLE_PRECISION)

warpVector = vtk.vtkWarpVector()
# warpVector.SetInputConnection(strainFilter.GetOutputPort())

# cell2point = vtk.vtkCellDataToPointData()
# cell2point.SetInputConnection(destroyTopology.GetOutputPort())
# cell2point.PassCellDataOff()
# cell2point.SetOutputPointsPrecision(vtk.vtkAlgorithm.DOUBLE_PRECISION)

plane = vtk.vtkPlane()
# plane.SetOrigin(0, 0, 0)
plane.SetNormal(0, 1, 0)

cutter = vtk.vtkCutter()
cutter.SetCutFunction(plane)
cutter.SetInputConnection(warpVector.GetOutputPort())
cutter.SetOutputPointsPrecision(vtk.vtkAlgorithm.DOUBLE_PRECISION)

writer = vtk.vtkXMLUnstructuredGridWriter()

ws = []
rs_contact = []
Fs = []

fig, ax = plt.subplots()

fig.subplots_adjust(right=0.75)
ax.set_xlabel(r"$\xi$ / m")
ax.set_ylabel(r"$\bar p$ / Pa")
add_leg = True


for t, fn in zip(ts, fns, strict=True):
    print("###### time", t)
    reader.SetFileName(fn)
    reader.Update()
    grid = reader.GetOutput()

    disp_2d = vtk_to_numpy(grid.GetPointData().GetArray("displacement"))
    disp_3d = np.zeros((disp_2d.shape[0], 3))
    disp_3d[:, (0, 1)] = disp_2d
    disp_3d_vtk = numpy_to_vtk(disp_3d, 1)
    disp_3d_vtk.SetName("u")

    grid.GetPointData().AddArray(disp_3d_vtk)
    grid.GetPointData().SetActiveVectors("u")
    warpVector.SetInputData(grid)
    warpVector.Update()
    grid = warpVector.GetOutput()

    xmin, xmax, ymin, ymax, zmin, zmax = grid.GetBounds()
    y_top = get_y_top(t)
    try:
        assert abs(ymax - y_top) < 1e-7
    except:
        print(ymax, y_top, ymax - y_top)
        raise
    ws.append(2 * (1.0 - y_top))

    # determine top boundary
    coords = vtk_to_numpy(grid.GetPoints().GetData())
    assert abs(min(coords[:, 0])) < 1e-8
    idcs_top_boundary = np.where(coords[:, 1] > y_top - 1e-7)[0]
    # print(idcs_top_boundary)
    assert len(idcs_top_boundary) != 0
    xs_top_boundary = coords[idcs_top_boundary, 0]
    idx_max = np.argmax(xs_top_boundary)
    r_contact = max(xs_top_boundary[idx_max], 0.0)
    y_at_r_contact = coords[idcs_top_boundary[idx_max], 1]
    print("radius of contact area:", r_contact, "at y =", y_at_r_contact)

    rs_contact.append(r_contact)

    def average_stress(rs, stress):
        r_contact = max(rs)
        rs_int = np.linspace(min(rs), max(rs) - 1e-8, max(len(rs), 200))
        stress_int = interp1d(rs, stress, bounds_error=False, fill_value=0.0)
        avg_stress = np.empty_like(rs_int)

        for i, r in enumerate(rs_int):
            rho_max = np.sqrt(r_contact**2 - r**2)
            rhos = np.linspace(0, rho_max, 100)
            xis = np.sqrt(rhos**2 + r**2)
            try:
                assert max(xis) <= r_contact + 1e-8
            except:
                print(max(xis), r_contact)
                raise
            avg_stress[i] = 1.0 / rho_max * np.trapezoid(x=rhos, y=stress_int(xis))
        # avg_stress[-1] = 0.0

        return rs_int, avg_stress

    def total_force(rs, stress):
        assert all(rs > -1e-6)
        rs_int = np.linspace(min(rs), max(rs), max(len(rs), 200))
        stress_int = interp1d(rs, stress, bounds_error=False, fill_value=0.0)

        return 2.0 * np.pi * np.trapezoid(x=rs_int, y=rs_int * stress_int(rs_int))

    if t > 0:
        plane.SetOrigin(0, y_at_r_contact, 0)
        cutter.Update()
        grid = cutter.GetOutput()

        # for a in range(grid.GetPointData().GetNumberOfArrays()):
        #     print(grid.GetPointData().GetArrayName(a))

        xs = vtk_to_numpy(grid.GetPoints().GetData())[:, 0]
        try:
            assert abs(max(xs) - r_contact) < 1e-7
            assert abs(min(xs)) < 1e-8
        except:
            print(min(xs), max(xs), r_contact, max(xs) - r_contact)
            raise

        w_0 = 2.0 * (1.0 - y_top)

        rs = np.linspace(0, r_contact, 200)
        if add_leg:
            ax.plot([], [], color="k", ls=":", label="ref")
        (h,) = ax.plot(rs, -p_contact(rs, r_contact), ls=":")

        stress_yy = vtk_to_numpy(grid.GetPointData().GetArray("sigma"))[:, 1]
        rs, avg_stress_yy = average_stress(xs, stress_yy)
        if add_leg:
            ax.plot([], [], color="k", ls="-", label="ogs")
        ax.plot(
            rs,
            avg_stress_yy,
            color=h.get_color(),
            ls="-",
            label=rf"$w_0 = {w_0}$",
        )

        ax.scatter([r_contact], [0], color=h.get_color())
        ax.legend(loc="center left", bbox_to_anchor=(1.0, 0.5))

        fig.savefig("stress_at_contact.png")
        add_leg = False

        Fs.append(-total_force(xs, stress_yy))

fig.savefig("stress_at_contact.png")
plt.close(fig)


fig, ax = plt.subplots()
ax.scatter(ws, rs_contact, label="ogs")

ws_ref = np.linspace(0, max(ws), 200)
rs_ref = np.sqrt(ws_ref * R)

ax.plot(ws_ref, rs_ref, label="ref")

ax.legend()

ax.set_xlabel("displacement of sphere centers $w_0$ / m")
ax.set_ylabel("radius of contact area $a$ / m")

fig.subplots_adjust(left=0.15, right=0.98)
fig.savefig("contact_radii.png")
plt.close(fig)


fig, ax = plt.subplots()
ax.scatter(rs_contact[1:], Fs, label="ogs")

rs_ref = np.linspace(0, max(rs_contact), 200)
Fs_ref = 8.0 * rs_ref**3 * kappa / 3.0

ax.plot(rs_ref, Fs_ref, label="ref")

legend = ax.legend()
legend.get_frame().set_facecolor("white")

ax.set_xlabel("radius of contact area $a$ / m")
ax.set_ylabel("applied force $F$ / N")

fig.subplots_adjust(left=0.15, right=0.98)
fig.savefig("total_force.png")
plt.close(fig)
