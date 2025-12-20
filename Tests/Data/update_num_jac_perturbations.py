import os
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import ogstools as ot
import pyvista as pv


def parse_space_separated_floats(text: str) -> np.ndarray:
    """Parse multiline text into a 2D NumPy array of floats."""
    # Split by lines, then split each line into floats
    rows = []
    for line in text.strip().splitlines():
        if line.strip():
            #            rows.append([float(x) for x in line.split()])
            for x in line.split():
                rows.append(float(x))
    return np.array(rows)


num_non_displacement_components = {
    "STEADY_STATE_DIFFUSION": 1,
    "RICHARDS_FLOW": 1,
    "RICHARDS_MECHANICS": 1,
    "THERMO_RICHARDS_MECHANICS": 2,
    "HEAT_CONDUCTION": 1,
    "HEAT_TRANSPORT_BHE": 2,
    "TWOPHASE_FLOW_PP": 2,
    "THERMAL_TWOPHASE_WITH_PP": 4,
    "TH2M": 3,
    "WELLBORE_SIMULATOR": 3,
}


def get_perturbation_data(prj_path: str):
    tree = ET.parse(prj_path)
    root = tree.getroot()
    jacobian_assembler = root.find(".//jacobian_assembler")
    if jacobian_assembler is None:
        return None

    jacobian_assembler_type = root.find(".//jacobian_assembler/type").text.strip()
    print(jacobian_assembler_type)
    if (jacobian_assembler_type != "CentralDifferences") and (
        jacobian_assembler_type != "ForwardDifferences"
    ):
        return None

    process_type = root.find(".//processes/process/type").text.strip()
    print(f"Process type: {process_type}")

    print(f"Processing project file {prj_path}")

    nndc = num_non_displacement_components[process_type]
    eps_tag = root.find(".//jacobian_assembler/epsilons")
    if eps_tag is not None and eps_tag.text:
        msg = "File already has <epsilons> tag."
        print(msg)
        return None

    component_mag_tag = root.find(".//component_magnitudes")
    rel_eps_tag = root.find(".//relative_epsilons")
    if (component_mag_tag is None or not component_mag_tag.text) and (
        rel_eps_tag is None or not rel_eps_tag.text
    ):
        return {
            "data_specified": False,
            "relative_epsilons": np.array([1.0e-8] * nndc),
        }
        # raise RuntimeError("No <component_magnitudes> tag or empty text found.")

    component_magnitudes = parse_space_separated_floats(component_mag_tag.text)
    print("Original component magnitudes:")
    print(component_magnitudes)

    relative_epsilons = parse_space_separated_floats(rel_eps_tag.text)
    print("Original relative epsilons:")
    print(relative_epsilons)
    if len(relative_epsilons) != len(component_magnitudes):
        msg = (
            f"The number of component magnitudes {component_magnitudes} is not equal to "
            f"the number of relative epsilons {relative_epsilons}"
        )
        raise RuntimeError(msg)

    if len(relative_epsilons) == nndc and len(component_magnitudes) == nndc:
        # print(f"File {prj_path} has already been updated.")
        return {
            "data_specified": True,
            "relative_epsilons": relative_epsilons * component_magnitudes,
        }

    first_mesh_tag = root.find(".//meshes/mesh")
    if first_mesh_tag is None:
        first_mesh_tag = root.find(".//mesh")
        if first_mesh_tag is None:
            msg = "No <mesh> tag found."
            raise RuntimeError(msg)

    mesh_filename = first_mesh_tag.text.strip()
    print(f"Mesh file name: {mesh_filename}")
    mesh = pv.read(mesh_filename)

    max_cell_points = max(mesh.get_cell(i).n_points for i in range(mesh.n_cells))
    print(f"Maximum number of nodes per element: {max_cell_points}")
    mesh_dimension = max(mesh.get_cell(i).dimension for i in range(mesh.n_cells))
    print(f"Mesh dimension: {mesh_dimension}")

    epsilons = relative_epsilons * component_magnitudes

    old_num_perburtations = (nndc + mesh_dimension) * max_cell_points
    if old_num_perburtations > len(relative_epsilons):
        print(
            f"The expected number of relative epsilons is {old_num_perburtations}."
            f"However only {len(relative_epsilons)} relative epsilons value are given."
        )
        return {
            "data_specified": True,
            "relative_epsilons": np.array([epsilons[0]] * nndc),
        }
    if old_num_perburtations < len(relative_epsilons):
        print(
            f"The expected number of relative epsilons is {old_num_perburtations}."
            f"However {len(relative_epsilons)} relative epsilons value are given."
            f"The extra values will be ignored."
        )
        return {
            "data_specified": True,
            "relative_epsilons": np.array([epsilons[0]] * nndc),
        }

    return {
        "data_specified": True,
        "relative_epsilons": np.array(
            [epsilons[i * max_cell_points] for i in range(nndc)]
        ),
    }


def update_perturbation(prj_file: str):
    data = get_perturbation_data(prj_file)

    if data is None:
        # No change needed
        return

    print(data["relative_epsilons"])

    # out_file = os.path.splitext(prj_file)[0]
    # print(out_file)
    prj = ot.Project(input_file=prj_file, output_file=Path("", f"{prj_file}"))

    eps_str = " ".join(map(str, data["relative_epsilons"]))

    if data["data_specified"]:
        prj.remove_element(
            "./processes/process/jacobian_assembler/component_magnitudes"
        )
        prj.remove_element("./processes/process/jacobian_assembler/relative_epsilons")
    prj.add_element(
        parent_xpath="./processes/process/jacobian_assembler",
        tag="epsilons",
        text=eps_str,
    )

    prj.write_input()


def process_prj_files(root_dir="."):
    """
    Recursively search all sub-directories for .prj files,
    change directory to where each file is located,
    call update_perturbation(prj_file),
    then continue searching other files.
    """
    # Remember the starting directory so we can return after each file
    start_dir = Path.cwd()

    for dirpath, _, filenames in os.walk(root_dir):
        for filename in filenames:
            if filename.lower().endswith(".prj"):
                # Change working directory to where the file is located
                os.chdir(dirpath)
                print(f"Changed directory to: {dirpath}")

                try:
                    update_perturbation(filename)
                finally:
                    # Always return to the starting directory
                    os.chdir(start_dir)


def main():
    if "--help" in sys.argv:
        print(
            "This script updates the perturbations of numerical Jacobian to the lastest OGS."
        )
        print("Usage: update_num_jac_perturbations.py")
        return

    # Start searching from the current directory
    process_prj_files()


if __name__ == "__main__":
    main()
