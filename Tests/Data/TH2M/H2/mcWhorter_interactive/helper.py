# SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
# SPDX-License-Identifier: BSD-3-Clause

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
import vtk
from IPython.display import HTML, Javascript, display


class plot_with_error:
    def __init__(self, x_ref, y_ref, style, label_ref, xlabel, ylabel, dylabel):
        self.x_ref = x_ref
        self.y_ref = y_ref
        self.label_ref = label_ref
        self.style = style

        self.xlabel = xlabel
        self.ylabel = ylabel
        self.dylabel = dylabel

        self.data = []

    def append(self, x, y, style, label):
        self.data.append([x, y, style, label])

    def plot(self, xlim=None):
        _fig, (ax0, ax1) = plt.subplots(2, 1, sharex=True, figsize=(6, 6))

        ax0.plot(self.x_ref, self.y_ref, self.style, label=self.label_ref)

        for data in self.data:
            ax0.plot(data[0], data[1], data[2], label=data[3])

        ax0.set_ylabel(self.ylabel)
        ax0.legend()
        ax0.yaxis.tick_left()

        # use default parameter in rcParams, not calling tick_right()
        for data in self.data:
            x, y, style, label = data
            dy = y - np.interp(x, self.x_ref, self.y_ref)
            ax1.plot(x, dy, style, label=label)

        ax1.set_xlabel(self.xlabel)
        ax1.set_ylabel(self.dylabel)

        if not None:
            ax1.set_xlim(xlim)

        plt.show()


def create_1d_mesh(point_a, point_b, num_points, mesh_type="line"):
    points = np.linspace(point_a, point_b, num_points)
    vtk_points = vtk.vtkPoints()
    cells = vtk.vtkCellArray()

    def add_point(point):
        for i in range(vtk_points.GetNumberOfPoints()):
            if np.allclose(vtk_points.GetPoint(i), point):
                return i
        return vtk_points.InsertNextPoint(point)

    if mesh_type == "line":
        for point in points:
            add_point(point)
        for i in range(num_points - 1):
            line = vtk.vtkLine()
            line.GetPointIds().SetId(0, i)
            line.GetPointIds().SetId(1, i + 1)
            cells.InsertNextCell(line)
        cell_type = vtk.VTK_LINE
    elif mesh_type == "quad":
        d = np.linalg.norm(np.array(point_b) - np.array(point_a))
        offset = np.array([0, d / num_points, 0])

        all_points = np.vstack((points, points + offset))
        point_ids = [add_point(p) for p in all_points]

        for i in range(num_points - 1):
            quad = vtk.vtkQuad()
            quad.GetPointIds().SetId(0, point_ids[i])
            quad.GetPointIds().SetId(1, point_ids[i + 1])
            quad.GetPointIds().SetId(2, point_ids[i + 1 + num_points])
            quad.GetPointIds().SetId(3, point_ids[i + num_points])
            cells.InsertNextCell(quad)
        cell_type = vtk.VTK_QUAD
    else:
        msg = f"Invalid mesh_type '{mesh_type}'. Choose 'line' or 'quad'."
        raise ValueError(msg)

    mesh = vtk.vtkUnstructuredGrid()
    mesh.SetPoints(vtk_points)
    mesh.SetCells(cell_type, cells)

    mesh = pv.UnstructuredGrid(mesh)
    mesh.point_data["bulk_node_ids"] = np.arange(0, mesh.n_points, dtype=np.uint64)
    mesh.cell_data["bulk_element_ids"] = np.arange(0, mesh.n_cells, dtype=np.uint64)

    return mesh


def create_boundary_line_meshes(point_a, point_b, num_points):
    d = np.linalg.norm(np.array(point_b) - np.array(point_a))
    offset = np.array([0, d / num_points, 0])

    def create_line_mesh(start_point):
        points = vtk.vtkPoints()
        points.InsertNextPoint(start_point)
        points.InsertNextPoint(start_point + offset)

        line = vtk.vtkLine()
        line.GetPointIds().SetId(0, 0)
        line.GetPointIds().SetId(1, 1)

        cells = vtk.vtkCellArray()
        cells.InsertNextCell(line)

        mesh = vtk.vtkPolyData()
        mesh.SetPoints(points)
        mesh.SetLines(cells)
        return pv.UnstructuredGrid(mesh)

    return tuple(map(create_line_mesh, [point_a, point_b]))


def format_numbers(x):
    if x == 0:
        return "0.0"
    if 0.0001 <= abs(x) <= 999999:
        return f"{x:.4f}".rstrip("0").rstrip(".")
    return f"{x:.4e}"


def render_latex_table(pandas_dataframe, latex_column):
    """
    Renders a Pandas DataFrame table in a Jupyter Notebook,
    ensuring that LaTeX formulas are correctly displayed.

    :param pandas_dataframe: Pandas DataFrame
    :param latex_column: Name of column with LaTeX equations
    """

    pandas_table = pandas_dataframe.copy()

    # Create HTML for LaTeX column
    pandas_table[latex_column] = pandas_table[latex_column].apply(
        lambda x: f"<div style='text-align: center;'>$$ {x} $$</div>"  # MathJax syntax
    )

    # Convert Dataframe to HTML and use MathJax to display equations
    html_table = pandas_table.to_html(
        escape=False, index=False, float_format=format_numbers
    )  # deactivate Escape to render HTML

    display(
        HTML(
            """
        <script type="text/javascript" async
            src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.7/MathJax.js?config=TeX-MML-AM_CHTML">
        </script>
        """
            + html_table
        )
    )

    display(
        Javascript(
            """
        if (window.MathJax) {
            MathJax.Hub.Queue(["Typeset", MathJax.Hub]);
        }
    """
        )
    )
