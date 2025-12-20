# SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
# SPDX-License-Identifier: BSD-3-Clause

from pathlib import Path
from subprocess import run

import numpy as np
import ogstools as ot
import pyvista as pv


class SingleOGSModel:
    """A single OGS run model"""

    def __init__(
        self,
        model,
        mesh_path,
        out_dir,
        output_prefix="",
        n_fracture_p_ncs=3,
        method="LIE",
        model_type="default",
        mesh_size=None,
        tension_cutoff=False,
        fracture_model_type=None,
        materials=None,
        n_mpi=1,
    ):
        self.model = model
        self.method = method
        self.n_fracture_p_ncs = n_fracture_p_ncs
        self.model_type = model_type
        self.mesh_size = mesh_size
        self.tension_cutoff = tension_cutoff
        self.fracture_model_type = fracture_model_type
        self.output_prefix = output_prefix
        self.out_dir = out_dir
        self.mesh_path = Path(self.out_dir, mesh_path)
        self.materials = materials
        self.current_load_case = "default"
        self.n_mpi = n_mpi

        print(f"mesh path: {mesh_path}")
        print(f"[DEBUG] Current working directory: {Path.cwd()}")
        print(f"[DEBUG] output project file: {self.model.prjfile}")

    def replace_curve(
        self,
        name=None,
        value=None,
        coords=None,
        parametertype=None,
        valuetag="values",
        coordstag="coords",
    ):
        root = self.model._get_root()
        parameterpath = "./curves/curve"
        parameterpointer = self.model._get_parameter_pointer(root, name, parameterpath)
        self.model._set_type_value(
            parameterpointer, value, parametertype, valuetag=valuetag
        )
        self.model._set_type_value(
            parameterpointer, coords, parametertype, valuetag=coordstag
        )

    def apply_load_curves(self, curve_name, coords, values):
        values2string = " ".join(map(str, values))
        self.replace_curve(name=curve_name, value=values2string, coords=coords)

    def apply_all_loads(self, times_string, pee_load_values):
        dss_load_values = pee_load_values.copy()
        ncols = dss_load_values.shape[1]
        for i in range(ncols):
            j = (i + 1) % ncols
            dss_load_values[:, i] = 0.5 * (
                dss_load_values[:, i] + dss_load_values[:, j]
            )

        for i in range(ncols):
            idx = i + 1
            pee_tag = f"PEE{idx}_SURFACE_CURVE"
            dss_tag = f"DSS{idx}_SURFACE_CURVE"

            self.apply_load_curves(
                pee_tag, times_string, pee_load_values[:, i].tolist()
            )
            self.apply_load_curves(
                dss_tag, times_string, dss_load_values[:, i].tolist()
            )

        for i in reversed(range(ncols)):
            idx = i + 1
            pee_tag_a = f"PEE{idx}a_SURFACE_CURVE"
            dss_tag_a = f"DSS{idx}a_SURFACE_CURVE"

            self.apply_load_curves(
                pee_tag_a, times_string, pee_load_values[:, i].tolist()
            )
            self.apply_load_curves(
                dss_tag_a, times_string, dss_load_values[:, i].tolist()
            )

    def run(self, times_string, PEE_load_values, material_name):
        try:
            self.apply_all_loads(times_string, PEE_load_values)
            if self.method == "LIE":
                self.apply_materials_for_LIE(material_name, model_type=self.model_type)
            elif self.method == "VPF":
                self.apply_materials_for_VPF(material_name, model_type=self.model_type)

            output_file_name_prefix = self.output_prefix + "_" + material_name
            self.model.replace_text(
                output_file_name_prefix, xpath="./time_loop/output/prefix"
            )
            self.model.write_input()

            log_file_name = Path(self.out_dir, f"log_{output_file_name_prefix}.txt")
            if self.method == "VPF":
                wrapper = f"mpirun --bind-to none -n {self.n_mpi}"
            else:
                wrapper = None

            self.model.run_model(
                logfile=log_file_name,
                wrapper=wrapper,
                args=f"-o {self.out_dir} -m {self.mesh_path}",
            )

            return Path(self.out_dir, f"{output_file_name_prefix}.pvd")
        except Exception as e:
            print(f"An error occurred run: {e}")

    def change_model_dt0_dt_min(self, dt0, dt_min, model_type):
        if model_type == "HM2a":
            self.model.replace_text(
                dt0,
                xpath="./time_loop/processes/process/time_stepping/initial_dt",
            )
            self.model.replace_text(
                dt_min,
                xpath="./time_loop/processes/process/time_stepping/minimum_dt",
            )

    def change_model_for_impermeable_sample(self, n_bcs):
        self.model.replace_text(
            "true", xpath="./processes/process/deactivate_matrix_in_flow"
        )
        for i in range(n_bcs + self.n_fracture_p_ncs, self.n_fracture_p_ncs, -1):
            xpath = f'./process_variables/process_variable[name="pressure"]/boundary_conditions/boundary_condition[{i}]'
            self.model.remove_element(xpath)

    def update_fracture_effective_stress(self, model_type, load_case):
        if model_type == "HM2b":
            value = "0.0 -2.0e6"
        elif model_type == "HM2a":
            if load_case.upper() in ["C"]:
                value = "0.0 -5.5e6"
            elif load_case.upper() in ["A", "B", "E", "F"]:
                value = "0.0 -4.0e6"
            else:
                msg = f"Invalid load case '{load_case}' for model type HM2a."
                raise ValueError(msg)
        else:
            return
        xpath = ".//parameters/parameter[name='fracture_effective_stress0']/values"
        try:
            self.model.replace_text(value, xpath=xpath)
        except Exception as e:
            msg = f"Failed to update 'fracture_effective_stress0': {e}"
            raise RuntimeError(msg) from e

    def run_for_all_samples(
        self,
        times_string,
        PEE_load_values,
        material_names,
        load_case="A",
        use_b_bar_value="false",
    ):
        vtu_file_names = []
        separate_line = "-" * 60
        for material_name in material_names:
            print(f"\n{separate_line}")
            print(f"* Running the simulation for sample {material_name}:")
            print(f"{separate_line}\n")

            if load_case in ["A", "B", "C", "E", "F"]:
                self.change_model_dt0_dt_min(0.001, 0.001, model_type=self.model_type)
                self.update_fracture_effective_stress(
                    model_type=self.model_type, load_case=load_case
                )

            self.model.replace_text(
                use_b_bar_value, xpath="./processes/process/use_b_bar"
            )
            pvd_file_name = self.run(times_string, PEE_load_values, material_name)

            if pvd_file_name is None:
                msg = f"Simulation failed for {material_name}, no PVD file generated."
                raise ValueError(msg)
            vtu_file_names.append(pvd_file_name)

        return vtu_file_names

    def apply_materials_for_LIE(self, material_name, model_type="default"):
        try:
            material = self.materials[material_name]
            rubber = material["rubber_sheath"]

            self.model.replace_parameter_value(
                name="E1", value=material["young_sample"]
            )
            self.model.replace_parameter_value(name="nu1", value=material["nu_sample"])
            self.model.replace_parameter_value(name="E2", value=rubber["young_modulus"])
            self.model.replace_parameter_value(
                name="nu2", value=rubber["poisson_ratio"]
            )

            if model_type in ["HM1", "HM2a", "HM2b"]:
                rate = material["fluid"]["injectionFlowRate_Inlet"]
                self.model.replace_parameter_value(
                    name="injectionFlowRate_Inlet",
                    value=rate,
                )

            if model_type in ["HM2a", "HM2b"]:
                pout = material["fluid"]["p_outlet"]
                self.model.replace_parameter_value(
                    name="p_outlet",
                    value=pout,
                )

            if model_type in ["M1", "M2a", "M2b"]:
                for mid in [0, 1]:
                    self.model.replace_phase_property_value(
                        mediumid=mid,
                        phase="Solid",
                        name="density",
                        value=material["density_solid"],
                    )
                self.model.replace_phase_property_value(
                    mediumid=2,
                    phase="Solid",
                    name="density",
                    value=rubber["density"],
                )

            if model_type in ["HM1", "HM2a", "HM2b"]:
                mids = [0, 1, 2] if model_type == "HM1" else [0, 1, 2, 3]
                for mid in mids:
                    self.model.replace_phase_property_value(
                        mediumid=mid,
                        phase="AqueousLiquid",
                        name="density",
                        value=material["fluid"]["density"],
                    )
                    self.model.replace_phase_property_value(
                        mediumid=mid,
                        phase="AqueousLiquid",
                        name="viscosity",
                        value=material["fluid"]["viscosity"],
                    )

                    if mid == 2:
                        # medium 2 is rubber
                        sd = (
                            rubber["density_solid"]
                            if "density_solid" in rubber
                            else rubber["density"]
                        )
                        src = "rubber"
                    else:
                        sd = (
                            material["density_solid"]
                            if "density_solid" in material
                            else material["density"]
                        )
                        src = "material"

                    self.model.replace_phase_property_value(
                        mediumid=mid,
                        phase="Solid",
                        name="density",
                        value=sd,
                    )
                    print(f"[mid={mid}] Solid density set from {src}: {sd}")

                    if mid != 3:
                        if mid < 2:
                            bc, perm, poro = (
                                material["biot"],
                                material["permeability"],
                                material["porosity"],
                            )
                        else:
                            bc, perm, poro = (
                                rubber["biot"],
                                rubber["permeability"],
                                rubber["porosity"],
                            )

                        self.model.replace_medium_property_value(
                            mediumid=mid,
                            name="biot_coefficient",
                            value=bc,
                        )
                        self.model.replace_medium_property_value(
                            mediumid=mid,
                            name="permeability",
                            value=perm,
                        )
                        self.model.replace_medium_property_value(
                            mediumid=mid,
                            name="porosity",
                            value=poro,
                        )
                        print(f"[mid={mid}] Bulk properties set")
                    else:
                        print(f"[mid={mid}], fracture, Bulk skipped")

            if self.n_fracture_p_ncs > 0:
                if model_type in ["HM2a", "HM2b"]:
                    self.model.replace_parameter_value(
                        name="S_f", value=material["S_f"]
                    )

                self.model.replace_parameter_value(name="Kn", value=material["k_n"])
                self.model.replace_parameter_value(name="Ks", value=material["k_t"])
                self.model.replace_parameter_value(
                    name="aperture0", value=material["w_init"]
                )

                if self.fracture_model_type == "CohesiveZoneModeI":
                    print("DEBUG: Applying CohesiveZoneModeI fracture model")
                    self.model.replace_text(
                        "CohesiveZoneModeI",
                        xpath="./processes/process/fracture_model/type",
                    )

                    if self.model.tree.xpath(
                        "./processes/process/fracture_model/fracture_toughness"
                    ):
                        self.model.replace_text(
                            "Gc",
                            xpath="./processes/process/fracture_model/fracture_toughness",
                        )
                    else:
                        self.model.add_element(
                            "./processes/process/fracture_model",
                            "fracture_toughness",
                            "Gc",
                        )

                    if self.model.tree.xpath(
                        "./processes/process/fracture_model/peak_normal_traction"
                    ):
                        self.model.replace_text(
                            "t_np",
                            xpath="./processes/process/fracture_model/peak_normal_traction",
                        )
                    else:
                        self.model.add_element(
                            "./processes/process/fracture_model",
                            "peak_normal_traction",
                            "t_np",
                        )

                    self.model.replace_parameter_value(name="Gc", value=material["Gc"])
                    self.model.replace_parameter_value(
                        name="t_np", value=material["t_np"]
                    )

                elif self.fracture_model_type == "LinearElasticIsotropic":
                    print("DEBUG: Applying LinearElasticIsotropic fracture model")
                    self.model.replace_text(
                        "LinearElasticIsotropic",
                        xpath="./processes/process/fracture_model/type",
                    )

                if not self.tension_cutoff:
                    self.model.replace_text(
                        0,
                        xpath="./processes/process/fracture_model/tension_cutoff",
                    )

        except KeyError as e:
            print(f"KeyError while applying LIE materials: {e}")
        except Exception as e:
            print(f"An error occurred while applying LIE materials: {e}")

    def apply_materials_for_VPF(self, material_name, model_type="default"):
        try:
            material = self.materials[material_name]
            rubber = material["rubber_sheath"]

            if self.fracture_model_type:
                self.model.replace_text(
                    self.fracture_model_type,
                    xpath="./processes/process/energy_split_model",
                )

            if model_type in ["M1", "M2b", "M2a"]:
                biot = 0.0
                permeability = 1e-10
                biot_rubber = 0.0
                permeability_rubber = 1e-10
            else:
                biot = material["biot"]
                permeability = material["permeability"]
                biot_rubber = rubber["biot"]
                permeability_rubber = rubber["permeability"]

            self.model.replace_parameter_value(
                name="E1", value=material["young_sample"]
            )
            self.model.replace_parameter_value(name="nu1", value=material["nu_sample"])
            self.model.replace_parameter_value(name="E2", value=rubber["young_modulus"])
            self.model.replace_parameter_value(
                name="nu2", value=rubber["poisson_ratio"]
            )

            if model_type in ["HM3d"]:
                self.model.replace_parameter_value(name="gc", value=material["Gc"])
                Q0 = material["fluid"]["injectionFlowRate_Inlet"]
                self.model.replace_parameter_value(
                    name="injectionFlowRate_Inlet",
                    value=str(Q0),
                )

            if model_type in ["HM2a", "HM2b"]:
                Q0 = material["fluid"]["injectionFlowRate_Inlet"] / self.mesh_size
                self.model.replace_parameter_value(
                    name="injectionFlowRate_Inlet", value=Q0
                )
                self.model.replace_parameter_value(
                    name="p_outlet", value=material["fluid"]["p_outlet"]
                )
            if model_type in ["HM1"]:
                Q0 = (
                    material["fluid"]["injectionFlowRate_Inlet"]
                    / material["fluid"]["density"]
                )
                self.model.replace_parameter_value(
                    name="injectionFlowRate_Inlet", value=Q0
                )

            for mid in [0, 1, 2]:
                self.model.replace_phase_property_value(
                    mediumid=mid,
                    phase="AqueousLiquid",
                    name="density",
                    value=material["fluid"]["density"],
                )
                self.model.replace_phase_property_value(
                    mediumid=mid,
                    phase="AqueousLiquid",
                    name="viscosity",
                    value=material["fluid"]["viscosity"],
                )

            for mid in [0, 1]:
                self.model.replace_medium_property_value(
                    mediumid=mid,
                    name="permeability",
                    value=permeability,
                )
            self.model.replace_medium_property_value(
                mediumid=2, name="permeability", value=permeability_rubber
            )

            if self.mesh_size:
                ls = 2 * self.mesh_size
                print(f"DEBUG: Setting ls to {ls} based on mesh size {self.mesh_size}")
                self.model.replace_parameter_value(name="ls", value=ls)
                print("ls", ls)

            self.model.replace_text(
                xpath="./processes/process/fluid_compressibility",
                value=material["c_f"],
            )

            for mid in [0, 1]:
                self.model.replace_phase_property_value(
                    mediumid=mid,
                    phase="Solid",
                    name="density",
                    value=material["density_solid"],
                )
                self.model.replace_phase_property_value(
                    mediumid=mid,
                    phase="Solid",
                    name="biot_coefficient",
                    value=biot,
                )
                self.model.replace_medium_property_value(
                    mediumid=mid, name="porosity", value=material["porosity"]
                )

            self.model.replace_phase_property_value(
                mediumid=2,
                phase="Solid",
                name="density",
                value=rubber["density"],
            )
            self.model.replace_medium_property_value(
                mediumid=2, name="porosity", value=rubber["porosity"]
            )
            self.model.replace_phase_property_value(
                mediumid=2,
                phase="Solid",
                name="biot_coefficient",
                value=biot_rubber,
            )

            self.model.write_input()
            print(
                f"Material properties for '{material_name}' (VPF) applied successfully."
            )

        except KeyError as e:
            print(f"KeyError while applying VPF materials: {e}")
        except Exception as e:
            print(f"An error occurred while applying VPF materials: {e}")

    @staticmethod
    def run_simulations_with_fracture(
        times,
        output_prefix,
        out_dir,
        base_project_file,
        mesh_path,
        load_cases,
        material_names,
        materials,
        n_fracture_p_ncs=3,
        method="LIE",
        crack_type="full",
        model_type="default",
        use_b_bar_value="false",
        mesh_size=0.0005,
        delta=0.00025,
        borehole_r=0.005,
        tension_cutoff=False,
        fracture_model_type=None,
        n_mpi=1,
    ):
        all_vtu_file_names = {key: [] for key in load_cases}

        if method == "VPF" and model_type not in ["M1", "HM1"]:
            mesh_full_path = Path(Path.cwd(), mesh_path)
            print(
                f"Modifying mesh for VPF method at {mesh_full_path} with crack type: {crack_type} and h1: {mesh_size}"
            )
            SingleOGSModel.set_crack_VPF(
                mesh_full_path,
                "domain.vtu",
                "domain.vtu",
                crack_type=crack_type,
                fracture_angle=0,
                h1=mesh_size,
                delta=delta,
                borehole_r=borehole_r,
                material_names=material_names,
                materials=materials,
                model_type=model_type,
                n_mpi=n_mpi,
            )

        separate_line = "=" * 60

        for key in load_cases:
            pee_load_values = np.array(
                [
                    np.zeros(len(load_cases[key]), dtype=float),
                    load_cases[key],
                    load_cases[key],
                ]
            )

            project = ot.Project(
                input_file=base_project_file,
                output_file=Path(out_dir, f"{output_prefix}_{key}.prj"),
            )

            sing_ogs_model = SingleOGSModel(
                model=project,
                mesh_path=mesh_path,
                out_dir=out_dir,
                output_prefix=output_prefix + "_" + key,
                n_fracture_p_ncs=n_fracture_p_ncs,
                method=method,
                model_type=model_type,
                mesh_size=mesh_size,
                tension_cutoff=tension_cutoff,
                fracture_model_type=fracture_model_type,
                materials=materials,
                n_mpi=n_mpi,
            )

            print(f"{separate_line}")
            print(
                f"Running simulation for load case: {key} with method: {method}, tension_cutoff: {tension_cutoff}"
            )
            print(f"{separate_line}")

            try:
                vtu_file_names = sing_ogs_model.run_for_all_samples(
                    times_string=times,
                    PEE_load_values=pee_load_values,
                    material_names=material_names,
                    load_case=key,
                    use_b_bar_value=use_b_bar_value,
                )
            except FileNotFoundError as e:
                print(f"Error: {e}")
                continue

            all_vtu_file_names[key] = vtu_file_names

        return all_vtu_file_names

    @staticmethod
    def set_crack_VPF(
        mesh_path,
        in_m,
        out_m,
        crack_type,
        fracture_angle,
        h1,
        delta,
        borehole_r,
        material_names,
        materials,
        model_type,
        n_mpi,
    ):
        chosen = material_names[0]
        w_init = materials[chosen]["w_init"]
        mesh = pv.read(Path(mesh_path, in_m))
        # --- pointwise data ---
        pts = mesh.points[:, :2]
        angle = np.radians(fracture_angle)
        R = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
        pts_rot = pts.dot(R.T)
        tol = 0.001 * h1

        mask_pts = np.zeros(len(pts), dtype=bool)
        if model_type == "HM3d":
            dist = np.linalg.norm(pts_rot, axis=1)
            mask_pts = dist <= (borehole_r + delta + h1 * tol)

        else:
            if crack_type == "full":
                mask_pts = (np.abs(pts_rot[:, 0]) <= (0.094 + tol)) & (
                    np.abs(pts_rot[:, 1]) < h1
                )
            elif crack_type == "half":
                mask_pts = (
                    (pts_rot[:, 0] >= -tol)
                    & (pts_rot[:, 0] <= (0.04 + tol))
                    & (np.abs(pts_rot[:, 1]) < h1)
                )

        phasefield = np.ones(len(pts))
        phasefield[mask_pts] = 0.0
        mesh.point_data["pf-ic"] = phasefield

        # --- cell-wise data ---
        n_cells = mesh.n_cells
        mat_ids = np.ones(n_cells, dtype=np.int32)
        width_ic = np.zeros(n_cells)
        centers = mesh.cell_centers().points[:, :2]
        pts_rot_c = centers.dot(R.T)
        x, y = centers[:, 0], centers[:, 1]

        mask_cells = np.zeros(n_cells, dtype=bool)
        if model_type == "HM3d":
            dist_c = np.linalg.norm(pts_rot_c, axis=1)
            mask_cells = dist_c <= (borehole_r + delta + h1 * tol)

        else:
            if crack_type == "full":
                mask_cells = (np.abs(pts_rot_c[:, 0]) <= (0.094 + tol)) & (
                    np.abs(pts_rot_c[:, 1]) < 0.5 * h1
                )
            elif crack_type == "half":
                mask_cells = (
                    (pts_rot_c[:, 0] >= -tol)
                    & (pts_rot_c[:, 0] <= (0.04 + tol))
                    & (np.abs(pts_rot_c[:, 1]) < 0.5 * h1)
                )

        width_ic[mask_cells] = w_init
        dist2 = x**2 + y**2
        r_inner = 0.065**2
        r_mid_sq = 0.097**2 - 0.001 * h1

        mask_mid = (dist2 > r_inner) & (dist2 < r_mid_sq)
        mask_outer = dist2 >= r_mid_sq

        mat_ids[mask_cells] = 0
        mat_ids[mask_mid] = 1
        mat_ids[mask_outer] = 2

        mesh.cell_data["MaterialIDs"] = mat_ids
        mesh.cell_data["width_ic"] = width_ic
        mesh.save(Path(mesh_path, out_m), binary=False)

        def run_partmesh(mesh_path, n_mpi, physical_groups):
            run(
                f"partmesh -s -o {mesh_path} -i {mesh_path}/domain.vtu",
                shell=True,
                check=True,
            )
            run(
                f"partmesh -m -n {n_mpi} -o {mesh_path} -i {mesh_path}/domain.vtu -- "
                + " ".join(f"{mesh_path}/{group}" for group in physical_groups),
                shell=True,
                check=True,
            )

        HM3D_GROUPS = [
            "physical_group_DSS1.vtu",
            "physical_group_DSS1a.vtu",
            "physical_group_DSS2.vtu",
            "physical_group_DSS2a.vtu",
            "physical_group_DSS3.vtu",
            "physical_group_DSS3a.vtu",
            "physical_group_DSS4.vtu",
            "physical_group_DSS4a.vtu",
            "physical_group_DSS5.vtu",
            "physical_group_DSS5a.vtu",
            "physical_group_DSS6.vtu",
            "physical_group_DSS6a.vtu",
            "physical_group_DSS7.vtu",
            "physical_group_DSS7a.vtu",
            "physical_group_DSS8.vtu",
            "physical_group_DSS8a.vtu",
            "physical_group_PEE1.vtu",
            "physical_group_PEE1a.vtu",
            "physical_group_PEE2.vtu",
            "physical_group_PEE2a.vtu",
            "physical_group_PEE3.vtu",
            "physical_group_PEE3a.vtu",
            "physical_group_PEE4.vtu",
            "physical_group_PEE4a.vtu",
            "physical_group_PEE5.vtu",
            "physical_group_PEE5a.vtu",
            "physical_group_PEE6.vtu",
            "physical_group_PEE6a.vtu",
            "physical_group_PEE7.vtu",
            "physical_group_PEE7a.vtu",
            "physical_group_PEE8.vtu",
            "physical_group_PEE8a.vtu",
            "physical_group_p_bottom.vtu",
            "physical_group_p_left.vtu",
            "physical_group_p_right.vtu",
            "physical_group_p_top.vtu",
            "physical_group_borehole_boundary.vtu",
        ]

        HM2_GROUPS = HM3D_GROUPS[:-1] + [  # same as HM3D but replace last element
            "physical_group_Inlet.vtu",
            "physical_group_Outlet_R_embeddedFracture.vtu",
            "physical_group_Outlet_R_fullFracture.vtu",
            "physical_group_Outlet_L_fullFracture.vtu",
        ]

        if model_type == "HM3d":
            run_partmesh(mesh_path, n_mpi, HM3D_GROUPS)
        elif model_type in ["M2a", "M2b", "HM2a", "HM2b"]:
            run_partmesh(mesh_path, n_mpi, HM2_GROUPS)
