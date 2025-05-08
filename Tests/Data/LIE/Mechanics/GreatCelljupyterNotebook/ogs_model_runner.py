from pathlib import Path

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

            log_file_name = Path(self.out_dir, f"log_{self.output_prefix}.txt")
            self.model.run_model(
                logfile=log_file_name, args=f"-o {self.out_dir} -m {self.mesh_path}"
            )

            return Path(self.out_dir, f"{output_file_name_prefix}.pvd")
        except Exception as e:
            print(f"An error occurred run: {e}")

    def change_model_dt0_dt_min(self, dt0, dt_min):
        self.model.replace_text(
            dt0, xpath="./time_loop/processes/process/time_stepping/initial_dt"
        )
        self.model.replace_text(
            dt_min, xpath="./time_loop/processes/process/time_stepping/minimum_dt"
        )

    def change_model_for_impermeable_sample(self, n_bcs):
        self.model.replace_text(
            "true", xpath="./processes/process/deactivate_matrix_in_flow"
        )
        for i in range(n_bcs + self.n_fracture_p_ncs, self.n_fracture_p_ncs, -1):
            xpath = f'./process_variables/process_variable[name="pressure"]/boundary_conditions/boundary_condition[{i}]'
            self.model.remove_element(xpath)

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
                self.change_model_dt0_dt_min(0.0001, 0.0001)

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

            if model_type in ["M1", "M2a", "M2b"]:
                for mid in [0, 1]:
                    self.model.replace_phase_property_value(
                        mediumid=mid,
                        phase="Solid",
                        name="density",
                        value=material["density_solid"],
                    )
                self.model.replace_phase_property_value(
                    mediumid=2, phase="Solid", name="density", value=rubber["density"]
                )

            if model_type == "HM1":
                for mid in [0, 1]:
                    self.model.replace_phase_property_value(
                        mediumid=mid,
                        phase="Solid",
                        name="density",
                        value=material["density_solid"],
                    )
                    self.model.replace_medium_property_value(
                        mediumid=mid, name="biot_coefficient", value=material["biot"]
                    )
                    self.model.replace_medium_property_value(
                        mediumid=mid,
                        name="permeability",
                        value=material["permeability"],
                    )
                    self.model.replace_medium_property_value(
                        mediumid=mid, name="porosity", value=material["porosity"]
                    )

                self.model.replace_phase_property_value(
                    mediumid=2, phase="Solid", name="density", value=rubber["density"]
                )
                self.model.replace_medium_property_value(
                    mediumid=2, name="permeability", value=rubber["permeability"]
                )
                self.model.replace_medium_property_value(
                    mediumid=2, name="porosity", value=rubber["porosity"]
                )
                self.model.replace_medium_property_value(
                    mediumid=2, name="biot_coefficient", value=rubber["biot"]
                )

            if self.n_fracture_p_ncs > 0:
                if model_type not in ["M2a", "M2b"]:
                    self.model.replace_parameter_value(
                        name="k_permeability_sample", value=material["permeability"]
                    )
                    self.model.replace_parameter_value(
                        name="k_permeability_rubber", value=rubber["permeability"]
                    )
                    self.model.replace_parameter_value(
                        name="phi", value=material["porosity"]
                    )
                    self.model.replace_parameter_value(
                        name="S_f", value=material["S_f"]
                    )

                self.model.replace_parameter_value(name="Kn", value=material["k_n"])
                self.model.replace_parameter_value(name="Ks", value=material["k_t"])

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
                        0, xpath="./processes/process/fracture_model/tension_cutoff"
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

            if model_type in ["M1", "M2a", "M2b", "HM1"]:
                for mid in [0, 1]:
                    self.model.replace_medium_property_value(
                        mediumid=mid,
                        name="permeability",
                        value=permeability,
                    )
                self.model.replace_medium_property_value(
                    mediumid=2, name="permeability", value=permeability_rubber
                )
            else:
                self.model.replace_parameter_value(
                    name="k_permeability_sample", value=permeability
                )
                self.model.replace_parameter_value(
                    name="k_permeability_rubber", value=permeability_rubber
                )

            if self.mesh_size:
                ls = 2 * self.mesh_size
                print(f"DEBUG: Setting ls to {ls} based on mesh size {self.mesh_size}")
                self.model.replace_parameter_value(name="ls", value=ls)
                print("ls", ls)

            self.model.replace_text(
                xpath="./processes/process/fluid_compressibility", value=material["c_f"]
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
                mediumid=2, phase="Solid", name="density", value=rubber["density"]
            )
            self.model.replace_medium_property_value(
                mediumid=2, name="porosity", value=rubber["porosity"]
            )
            self.model.replace_phase_property_value(
                mediumid=2, phase="Solid", name="biot_coefficient", value=biot_rubber
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
        tension_cutoff=False,
        fracture_model_type=None,
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
            )

        separate_line = "=" * 60

        for key in load_cases:
            pee_load_values = -1 * np.array(
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
    def set_crack_VPF(mesh_path, in_m, out_m, crack_type, fracture_angle, h1):
        m = pv.read(Path(mesh_path, in_m))
        n_points = m.number_of_points
        pt_coords = m.points
        MaterialIDs = m.cell_data["MaterialIDs"]
        phasefield = np.ones(n_points)
        pressure = np.ones(n_points)

        r1 = 0.094  # full crack
        r2 = 0.04  # half crack

        angle_rad = np.radians(fracture_angle)

        for point_id in range(n_points):
            x = pt_coords[point_id, 0]
            y = pt_coords[point_id, 1]

            # Rotate the coordinates
            x_rotated = x * np.cos(angle_rad) - y * np.sin(angle_rad)
            y_rotated = x * np.sin(angle_rad) + y * np.cos(angle_rad)

            if crack_type == "full":
                if (-r1 - 0.001 * h1 <= x_rotated <= r1 + 0.001 * h1) and (
                    np.abs(y_rotated) < h1
                ):
                    phasefield[point_id] = 0
                    pressure[point_id] = 1.0e6
                else:
                    phasefield[point_id] = 1
                    pressure[point_id] = 1.0e5
            elif crack_type == "half":
                if (0 - 0.001 * h1 < x_rotated < r2 + 0.001 * h1) and (
                    np.abs(y_rotated) < h1
                ):
                    phasefield[point_id] = 0
                    pressure[point_id] = 1.0e6
                else:
                    phasefield[point_id] = 1
                    pressure[point_id] = 1.0e5
            else:
                phasefield = np.ones(n_points)
                pressure[point_id] = 1.0e5

        n_cells = m.n_cells
        width_ic = np.zeros(n_cells)

        for j in range(n_cells):
            cell = m.GetCell(j)
            x_min, x_max, y_min, y_max, z_min, z_max = cell.GetBounds()
            x = (x_min + x_max) / 2
            y = (y_min + y_max) / 2

            if crack_type == "full":
                if (-r1 - 0.001 * h1 <= x <= r1 + 0.001 * h1) and (
                    np.abs(y) < 0.5 * h1
                ):
                    width_ic[j] = 1e-5
                else:
                    width_ic[j] = 0.0
            elif crack_type == "half":
                if (0 <= x <= r2) and (np.abs(y) < 0.5 * h1):
                    width_ic[j] = 1e-5
                else:
                    width_ic[j] = 0.0

            # Material ID
            if (x**2 + y**2) <= 0.065**2:
                MaterialIDs[j] = 0
            elif 0.065**2 < (x**2 + y**2) <= 0.097**2 - 1e-3 * h1:
                MaterialIDs[j] = 1
            else:
                MaterialIDs[j] = 2

        m.cell_data["width_ic"] = width_ic
        m.point_data["pf-ic"] = phasefield
        m.point_data["p-ic"] = pressure
        m.cell_data["MaterialIDs"] = MaterialIDs

        m.save(Path(mesh_path, out_m), binary=False)
