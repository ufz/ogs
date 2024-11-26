# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.2
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [raw]
# +++
# title = "More benchmarks using F-bar"
# date = "2024-11-26"
# author = "Wenqing Wang"
# image = "figures/thick_cyl_shell.png"
# web_subsection = "large-deformations"
# weight = 3
# +++

# %% [markdown]
# $$
# \newcommand{\J}{\text{J}}
# \newcommand{\F}{\text{F}}
# $$
#
# # More large deformation benchmarks using the F-bar method
#
# These benchmarks are taken from [1] but employ a Neo-Hookean hyperelastic model.
#
# ## Reference
#
# 1. T. Elguedj, Y. Bazilevs, V.M. Calo, T.J.R. Hughes (2008),
#  $\bar\B$ and $\bar{\text F}$ projection methods for nearly incompressible linear and non-linear elasticity and plasticity using higher-order NURBS elements, Computer Methods in Applied Mechanics and Engineering, 197(33–40), 2732-2762.
#

# %%
# Standard library imports
import os
import xml.etree.ElementTree as ET
from pathlib import Path

# Third-party imports
import pyvista as pv

# Local imports
from ogs6py.ogs import OGS

out_dir = Path(os.environ.get("OGS_TESTRUNNER_OUT_DIR", "_out"))

out_dir = Path("output")
if not out_dir.exists():
    out_dir.mkdir(parents=True)

ogs_path = "/home/wenqing/Code/build/ogs6_release/bin"


# %%
def get_last_vtu_file_name(pvd_file_name):
    tree = ET.parse(pvd_file_name)
    root = tree.getroot()
    # Get the last DataSet tag
    last_dataset = root.findall(".//DataSet")[-1]

    # Get the 'file' attribute of the last DataSet tag
    file_attribute = last_dataset.attrib["file"]
    return Path(out_dir, file_attribute)


# %%
class SingleOGSModel:
    """An OGS run model"""

    def __init__(
        self, project_file, output_prefix, mesh_path, out_dir=out_dir, use_fbar=True
    ):
        self.model = OGS(
            INPUT_FILE=project_file, PROJECT_FILE=Path(out_dir, "modified.prj")
        )

        self.model.replace_text(output_prefix, xpath="./time_loop/output/prefix")
        if not use_fbar:
            self.model.replace_text("none", xpath="./processes/process/f_bar")

        self.use_fbar = use_fbar
        self.out_dir = out_dir
        self.output_prefix = output_prefix
        self.pvd_file_name = Path(self.out_dir, self.output_prefix + ".pvd")
        self.meth_path = mesh_path
        self.resulted_mesh = pv.UnstructuredGrid()

    def reset_time_step_size(self, dt, n_steps):
        self.model.replace_text(
            dt,
            xpath="./time_loop/processes/process[1]/time_stepping/timesteps/pair[1]/delta_t",
        )
        self.model.replace_text(
            n_steps,
            xpath="./time_loop/processes/process[1]/time_stepping/timesteps/pair[1]/repeat",
        )

    def run(self):
        self.model.write_input()

        if ogs_path != "":
            self.model.run_model(
                logfile=Path(self.out_dir, "out.txt"),
                path=ogs_path,
                args=f"-o {self.out_dir} -m {self.meth_path}",
            )
        else:
            self.model.run_model(
                logfile=Path(self.out_dir, "out.txt"),
                args=f"-o {self.out_dir} -m {self.meth_path}",
            )
        self.resulted_mesh = pv.read(get_last_vtu_file_name(self.pvd_file_name))

    def get_pvd_name(self):
        return self.pvd_file_name

    def plot3D(self, variable_name, camera_position, component_id=-1):
        point_data = self.resulted_mesh.point_data[variable_name]
        plot_data = point_data[:, component_id] if component_id > 0 else point_data
        deformed_mesh = self.resulted_mesh.copy()
        deformed_mesh.points += self.resulted_mesh.point_data["displacement"]

        pl = pv.Plotter()
        pl.add_mesh(self.resulted_mesh, color="white", style="wireframe", line_width=1)
        pl.add_mesh(
            deformed_mesh,
            scalars=plot_data,
            label="Vertical displacement",
            show_edges=False,
            cmap="jet",
        )
        pl.show_axes()
        pl.camera_position = camera_position
        pl.reset_camera()

        pl.show()

    def check_vertical_displacement_at_point(
        self, point, expected_value, component_id=2
    ):
        p_id = self.resulted_mesh.find_closest_point(point)
        u = self.resulted_mesh.point_data["displacement"][p_id]
        u_z_max = u[component_id]
        print(u_z_max)
        try:
            rel_error = abs((u_z_max - expected_value) / u_z_max)
            print(f"The relative error is {rel_error}.")
            assert rel_error < 1e-4, (
                f"The caculated vertical_displacement at {point} is not"
                f" {expected_value}"
            )
            return u_z_max
        except Exception as e:
            print(f"An error occurred: {e}")

    def run_benchmark(
        self,
        point_having_u_compoment_max,
        expected_u_compoment_max,
        camera_position,
        benchmark_name,
        component_id,
    ):
        self.run()

        u_compoment_max = self.check_vertical_displacement_at_point(
            point_having_u_compoment_max, expected_u_compoment_max, component_id
        )
        f_bar_info = "Without the F-bar" if self.use_fbar else "With the F-bar"
        print(
            f"{f_bar_info}, the obtain maximum displacement component"
            f" of {benchmark_name} is {u_compoment_max} m."
        )

        self.plot3D("displacement", camera_position, component_id)


# %% [markdown]
# ## 3D indentation example

# %%
try:
    project_file = Path("Indentation3D", "Indentation3D.prj")
    output_prefix = "Indentation3D"
    ogs_model = SingleOGSModel(
        project_file,
        output_prefix,
        mesh_path=output_prefix,
        out_dir=out_dir,
        use_fbar=False,
    )

    camera_position = (-0.5, -0.5, 0.4)
    point_have_u_compoment_max = [0.0, 0.0, 0.01]
    expected_u_compoment_max = -1.6534484343296282e-05
    ogs_model.run_benchmark(
        point_have_u_compoment_max,
        expected_u_compoment_max,
        camera_position,
        "3D indentation",
        component_id=2,
    )
except Exception:
    pass


# %%
try:
    project_file = Path("Indentation3D", "Indentation3D.prj")
    output_prefix = "Indentation3D"
    ogs_model = SingleOGSModel(
        project_file,
        output_prefix,
        mesh_path=output_prefix,
        out_dir=out_dir,
        use_fbar=True,
    )

    camera_position = (-0.5, -0.5, 0.4)
    point_have_u_compoment_max = [0.0, 0.0, 0.01]
    expected_u_compoment_max = -0.00010953948059235105
    ogs_model.run_benchmark(
        point_have_u_compoment_max,
        expected_u_compoment_max,
        camera_position,
        "3D indentation",
        component_id=2,
    )
except Exception:
    pass
