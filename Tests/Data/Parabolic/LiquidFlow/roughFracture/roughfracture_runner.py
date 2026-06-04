from __future__ import annotations

import shutil
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import ogstools as ot
import pyvista as pv


@dataclass
class RoughFractureCase:
    jrc: int
    sigma_mpa: float
    case_dir: Path
    mesh_dir: Path
    results_dir: Path
    mesh_file: Path
    prj_file: Path
    output_prefix: str
    ok: bool = False
    error: str | None = None


def sigma_tag(sigma: float) -> str:
    return str(float(sigma)).replace(".", "p")


def jrc_dir_name(jrc: int | float) -> str:
    return f"JRC_{int(float(jrc))}"


def print_mesh_inventory(src_mesh_dir: str | Path) -> None:
    src_mesh_dir = Path(src_mesh_dir)
    print(f"Input mesh directory  : {src_mesh_dir}")
    vtus = sorted(src_mesh_dir.glob("*.vtu"))
    print(f"Input mesh files      : {len(vtus)}")
    for vtu in vtus:
        print(f"  {vtu.name}")


def _update_density_property(project: ot.Project, density_spec: dict) -> None:
    if density_spec["type"] != "Linear":
        msg = (
            "The density property in roughFracture_synthesis_LF.prj is not "
            "Linear; reference_value, reference_condition, and slope are "
            "sub-elements of the Linear type only, so patching them for any "
            "other type silently does nothing."
        )
        raise ValueError(msg)

    prop_xpath = ".//media/medium/phases/phase[type='AqueousLiquid']/properties/property[name='density']"
    project.replace_phase_property_value(
        mediumid=0,
        phase="AqueousLiquid",
        name="density",
        value=density_spec["reference_value"],
        propertytype="Linear",
        valuetag="reference_value",
    )
    project.replace_text(
        str(density_spec["variable_name"]),
        xpath=f"{prop_xpath}/independent_variable/variable_name",
    )
    project.replace_text(
        str(density_spec["reference_condition"]),
        xpath=f"{prop_xpath}/independent_variable/reference_condition",
    )
    project.replace_text(
        str(density_spec["slope"]),
        xpath=f"{prop_xpath}/independent_variable/slope",
    )


def update_project_parameters(project: ot.Project, params: dict) -> None:
    project.set(output_prefix=params["prefix"], t_end=params["t_end"])
    project.replace_text(params["prefix"], xpath="./time_loop/output/prefix")
    project.replace_text(
        params["t_end"], xpath="./time_loop/processes/process/time_stepping/t_end"
    )
    project.replace_text(
        params["initial_dt"],
        xpath="./time_loop/processes/process/time_stepping/initial_dt",
    )
    project.replace_text(
        params["minimum_dt"],
        xpath="./time_loop/processes/process/time_stepping/minimum_dt",
    )
    project.replace_text(
        params["maximum_dt"],
        xpath="./time_loop/processes/process/time_stepping/maximum_dt",
    )
    project.replace_text(
        params["specific_body_force"], xpath="./processes/process/specific_body_force"
    )
    project.replace_text(
        params["equation_balance_type"],
        xpath="./processes/process/equation_balance_type",
    )
    project.replace_parameter_value("p0", params["initial_pressure"])
    project.replace_parameter_value("p_right", params["outlet_pressure"])
    project.replace_parameter_value("p_left", params["inlet_pressure"])
    project.replace_parameter_value(
        "constant_porosity_parameter", params["porosity_value"]
    )
    project.replace_parameter(
        "kappa1_frac",
        "MeshElement",
        ["field_name"],
        [params["permeability_mesh_field"]],
    )
    project.replace_parameter(
        "fracture_thickness_const",
        "MeshElement",
        ["field_name"],
        [params["fracture_thickness_mesh_field"]],
    )
    _update_density_property(project, params["fluid_density"])


def _ensure_aperture_and_cubic_perm(
    mesh_path: Path,
    *,
    ap_name: str = "aperture_closed",
    perm_name: str = "permeability_cubic",
    w_min: float = 1e-8,
) -> None:
    mesh = pv.read(mesh_path)
    if ap_name not in mesh.cell_data:
        msg = f"Cell field '{ap_name}' not found in {mesh_path}"
        raise KeyError(msg)
    aperture = np.maximum(np.asarray(mesh.cell_data[ap_name], float), float(w_min))
    mesh.cell_data[ap_name] = aperture
    mesh.cell_data[perm_name] = aperture**2 / 12.0
    if "MaterialIDs" not in mesh.cell_data:
        mesh.cell_data["MaterialIDs"] = np.zeros(mesh.n_cells, dtype=np.int32)
    mesh.field_data.clear()
    mesh.save(mesh_path)


def _extract_and_split_boundaries_2d(mesh_dir: Path, mesh_name_vtu: str) -> None:
    mesh = pv.read(mesh_dir / mesh_name_vtu)
    meshes = ot.Meshes.from_mesh(
        mesh,
        threshold_angle=None,
        domain_name=Path(mesh_name_vtu).stem,
    )
    meshes.save(mesh_dir, overwrite=True)


def prepare_case(
    *,
    jrc: int,
    sigma_mpa: float,
    src_mesh_dir: Path,
    prj_template: Path,
    base_run_dir: Path,
    user_parameters: dict,
    w_min: float = 1e-8,
) -> RoughFractureCase:
    s_tag = sigma_tag(sigma_mpa)
    src_mesh = src_mesh_dir / f"joint_JRC{int(jrc)}_sigma_{s_tag}MPa.vtu"
    if not src_mesh.exists():
        raise FileNotFoundError(src_mesh)

    case_dir = base_run_dir / jrc_dir_name(jrc) / f"sigma_{s_tag}"
    mesh_dir = case_dir / "mesh"
    results_dir = case_dir / "out"
    mesh_dir.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)

    mesh_file = mesh_dir / src_mesh.name
    shutil.copy2(src_mesh, mesh_file)
    _ensure_aperture_and_cubic_perm(mesh_file, w_min=w_min)
    _extract_and_split_boundaries_2d(mesh_dir, mesh_file.name)

    params = dict(user_parameters)
    output_prefix = f"{params['prefix']}_JRC{int(jrc)}_sigma_{s_tag}"
    params["prefix"] = output_prefix
    params["permeability_mesh_field"] = "permeability_cubic"
    params["fracture_thickness_mesh_field"] = "aperture_closed"

    prj_file = results_dir / prj_template.name.replace(".prj", "_final.prj")
    project = ot.Project(input_file=prj_template, output_file=prj_file)
    project.replace_mesh("joint_sigma_20MPa_midplane_tri.vtu", mesh_file.name)
    update_project_parameters(project, params)
    project.write_input()

    return RoughFractureCase(
        jrc=int(jrc),
        sigma_mpa=float(sigma_mpa),
        case_dir=case_dir,
        mesh_dir=mesh_dir,
        results_dir=results_dir,
        mesh_file=mesh_file,
        prj_file=prj_file,
        output_prefix=output_prefix,
    )


def run_case(case: RoughFractureCase) -> RoughFractureCase:
    try:
        project = ot.Project(input_file=case.prj_file, output_file=case.prj_file)
        project.run_model(
            args=f"-o {case.results_dir} -m {case.mesh_dir}",
            logfile=case.results_dir / "run.log",
        )
        if project.process.returncode != 0:
            msg = (
                f"OGS failed with return code {project.process.returncode}; "
                f"see {case.results_dir / 'run.log'}"
            )
            raise RuntimeError(msg)
        pvd_file = case.results_dir / f"{case.output_prefix}.pvd"
        if not pvd_file.exists():
            msg = f"Expected OGS output file was not created: {pvd_file}"
            raise FileNotFoundError(msg)
        case.ok = True
        case.error = None
    except Exception as exc:
        case.ok = False
        case.error = f"{type(exc).__name__}: {exc}"
    return case


def prepare_all_cases(
    *,
    out_dir: Path,
    base_prj: Path,
    src_mesh_dir: Path,
    user_parameters: dict,
    jrc_list: list[int],
    sigmas_mpa: list[float],
    w_min: float = 1e-8,
) -> list[RoughFractureCase]:
    base_run_dir = Path(out_dir) / "runs"
    base_run_dir.mkdir(parents=True, exist_ok=True)
    cases = []
    for jrc in jrc_list:
        for sigma in sigmas_mpa:
            case = prepare_case(
                jrc=int(jrc),
                sigma_mpa=float(sigma),
                src_mesh_dir=Path(src_mesh_dir),
                prj_template=Path(base_prj),
                base_run_dir=base_run_dir,
                user_parameters=user_parameters,
                w_min=w_min,
            )
            print(
                f"  prepared  JRC={case.jrc}, sigma={case.sigma_mpa:g} MPa -> {case.prj_file}"
            )
            cases.append(case)
    return cases


def run_all_cases(cases: list[RoughFractureCase]) -> list[RoughFractureCase]:
    results = []
    for case in cases:
        print(f"  running  JRC={case.jrc}, sigma={case.sigma_mpa:g} MPa")
        result = run_case(case)
        print(f"    {'OK' if result.ok else 'FAILED: ' + str(result.error)}")
        results.append(result)
    failed = [result for result in results if not result.ok]
    if failed:
        messages = "\n".join(
            f"JRC={case.jrc}, sigma={case.sigma_mpa:g} MPa: {case.error}"
            for case in failed
        )
        msg = f"{len(failed)} rough fracture OGS run(s) failed:\n{messages}"
        raise RuntimeError(msg)
    return results
