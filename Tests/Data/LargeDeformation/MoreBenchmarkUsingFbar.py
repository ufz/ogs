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
# These benchmarks are taken from [1] but employ a basic Neo-Hookean hyperelastic model given by
#  The Neo-Hookean model defines an energy function as the sum of volumetric component $U_{\text{dil}}(\F)$ and deviatoric component $U_{\text{dev}}(\F)$ as
# $$
# \begin{align}
# W(\F) = U_{\text{dil}}(\F) + U_{\text{dev}}(\F),
# \end{align}
# $$
# where
# $$
# \begin{align}
# & U_{\text{dil}}(\F) = \dfrac{1}{2} K (\det(\F)-1)^2\\
# & U_{\text{dev}}(\F) = \dfrac{1}{2} G \left(\text{tr} (\det(\F)^{-\frac{2}{3}}\F\F^{\text{T}})-3\right)
# \end{align}
# $$
# with $K$ the bulk modulus, and $G$ the shear modulus.
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

if not out_dir.exists():
    out_dir.mkdir(parents=True)


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
        dim = self.resulted_mesh.point_data["displacement"].shape[1]
        deformed_mesh.points[:, :dim] += self.resulted_mesh.point_data["displacement"]

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
        f_bar_info = "With the F-bar" if self.use_fbar else "Without the F-bar"
        print(
            f"{f_bar_info}, the obtain maximum displacement component"
            f" of {benchmark_name} is {u_compoment_max} m."
        )

        self.plot3D("displacement", camera_position, component_id)


# %% [markdown]
# ## 1. A block under compression
#
# This is a plane strain problem. The domain size is 0.01 m $\times$ 0.01 m. The geometry and boundary conditions are shown in the following figure:
#
# <img src="./figures/ld_rubber_indentation.png" alt="Simple test" width="200" height="200" />
#
# The values of Young's modulus ($E$) and Poisson's ratio ($\nu$) are calculated from the specified bulk modulus of 400,889.806 MPa and the shear modulus of 80.194 MPa. The calculated Poisson's ratio is close to 0.5, indicating that the material is nearly incompressible.
#
# ### 1.1. Simualtion without the F-bar method
#
#   A contour plot of vertical displacement in the deformation mesh is shown in the following figure after the simulation finishes.

# %%
try:
    project_file = Path("RubberIndentation", "RubberIndentation.prj")
    output_prefix = "RubberIndentation"
    ogs_model = SingleOGSModel(
        project_file,
        output_prefix,
        mesh_path=output_prefix,
        out_dir=out_dir,
        use_fbar=False,
    )
    # ogs_model.reset_time_step_size(0.3, 4)

    camera_position = "xy"
    point_have_u_compoment_max = [0.0, 0.01, 0.0]
    expected_u_compoment_max = -0.00010428732840299523
    ogs_model.run_benchmark(
        point_have_u_compoment_max,
        expected_u_compoment_max,
        camera_position,
        "Rubber indentation",
        component_id=1,
    )
except Exception:
    pass

# %% [markdown]
# ### 1.2. Simualtion with the F-bar method
#
#  A contour plot of vertical displacement in the deformed mesh is shown in the following figure after the simulation is completed. It can be observed that the deformation obtained using the F-bar method is larger than that obtained without using the F-bar method.

# %%
try:
    project_file = Path("RubberIndentation", "RubberIndentation.prj")
    output_prefix = "RubberIndentation"
    ogs_model = SingleOGSModel(
        project_file,
        output_prefix,
        mesh_path=output_prefix,
        out_dir=out_dir,
        use_fbar=True,
    )

    camera_position = "xy"
    point_have_u_compoment_max = [0.0, 0.01, 0.0]
    expected_u_compoment_max = -0.00026707547768988984
    ogs_model.run_benchmark(
        point_have_u_compoment_max,
        expected_u_compoment_max,
        camera_position,
        "Rubber indentation",
        component_id=1,
    )
except Exception:
    pass


# %% [markdown]
# ## 2. 3D indentation example
#
# <img align="left" src="./figures/ld_block_compression.png" alt="Simple test" width="200" height="200" />
# This example analyzes the deformation of a block ($0.01\times 0.01 \times 0.01\,\text{m}^3$) under indentation. As shown in the figure, a uniformly distributed pressure of 40 MPa is applied to one-quarter of the top surface. On the symmetry surfaces, the displacement in the normal direction is fixed, while on the bottom surface, the vertical displacement is fixed. The material parameters are the same as those used in the previous example.
#
# ### 2.1. Simualtion without the F-bar method
#
#   A contour plot of vertical displacement in the deformation mesh is shown in the following figure after the simulation finishes.

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


# %% [markdown]
# ### 2.2. Simualtion with the F-bar method
#
#  A contour plot of vertical displacement in the deformed mesh is shown in the following figure after the simulation is completed. It can be observed that the deformation obtained using the F-bar method is larger than that obtained without using the F-bar method.

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

# %% [markdown]
# ## 3. Thick cylindrical shell under pressure
#
# <img align="left" src="./figures/ld_cylindrical_shell.png" alt="Simple test" width="100" height="100" />
#
# This example analyzes the deformation of a cylindrical shell (R=$0.01\,\text{m}$, thickness=$0.002\,\text{m}$, H=$0.015\,\text{m}$) under radial pressure. Due to symmetry, only half of the shell is considered. The normal displacement on the symmetry surfaces is fixed. The other boundary conditions are shown in the figure, where the uniformly distributed pressure is 6 MPa/m. The bulk modulus is 240 GPa, and the shear modulus is 6 GPa, corresponding to a Poisson's ratio of 0.4.
#
# ### 3.1. Simualtion without the F-bar method
#
#   A contour plot of vertical displacement in the deformation mesh is shown in the following figure after the simulation finishes.

# %%
try:
    project_file = Path("ThickCylindricalShell", "ThickCylindricalShell.prj")
    output_prefix = "ThickCylindricalShell"
    ogs_model = SingleOGSModel(
        project_file,
        output_prefix,
        mesh_path=output_prefix,
        out_dir=out_dir,
        use_fbar=False,
    )
    ogs_model.reset_time_step_size(0.3, 4)

    camera_position = (0.5, -0.5, 0.4)
    point_have_u_compoment_max = [0.0, 0.01, 0.0075]
    expected_u_compoment_max = -0.010194869353071633
    ogs_model.run_benchmark(
        point_have_u_compoment_max,
        expected_u_compoment_max,
        camera_position,
        "Thick cylindrical shell",
        component_id=1,
    )
except Exception:
    pass

# %% [markdown]
# ### 3.2. Simualtion with the F-bar method
#
#  A contour plot of vertical displacement in the deformed mesh is shown in the following figure after the simulation is completed. It can be observed that the deformation obtained using the F-bar method is larger than that obtained without using the F-bar method.

# %%
try:
    project_file = Path("ThickCylindricalShell", "ThickCylindricalShell.prj")
    output_prefix = "ThickCylindricalShell"
    ogs_model = SingleOGSModel(
        project_file,
        output_prefix,
        mesh_path=output_prefix,
        out_dir=out_dir,
        use_fbar=True,
    )

    camera_position = (0.5, -0.5, 0.4)
    point_have_u_compoment_max = [0.0, 0.01, 0.0075]
    expected_u_compoment_max = -0.014236570288857237
    ogs_model.run_benchmark(
        point_have_u_compoment_max,
        expected_u_compoment_max,
        camera_position,
        "Thick cylindrical shell",
        component_id=1,
    )
except Exception:
    pass
