import colorsys
import math
from pathlib import Path
from typing import ClassVar

import matplotlib.pyplot as plt
import numpy as np
import ogstools as ot


class Plotter:
    DEFAULT_RC_PARAMS: ClassVar[dict] = {
        "font.family": "serif",
        "font.serif": ["DejaVu Serif"],
        "font.size": 12,
        "axes.labelsize": 12,
        "axes.titlesize": 14,
        "legend.fontsize": 8,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "lines.linewidth": 1.5,
        "axes.linewidth": 0.8,
        "grid.linewidth": 0.4,
        "savefig.dpi": 600,
        "mathtext.fontset": "stix",
    }

    METHOD_COLORS: ClassVar[dict] = {
        "VPF": "#27ae60",
        "LBLN-VPF-MOOSES": "#2980b9",
        "LIE": "#c0392b",
        "HM": "#c0392b",
        "SD": "#c0392b",
        "SNU-ML-PINN": "#d35400",
        "KIGAM": "#1a5276",
        "KIGAM-OGS3DEC": "#1f618d",
        "KAERI-TOUGH-3DEC": "#e67e22",
        "FEM": "#8e44ad",
        "Unknown": "#7f8c8d",
        "Default": "#7f8c8d",
    }

    MATERIAL_SHADES: ClassVar[dict] = {
        "Gneiss": 0.9,
        "Granite": 0.8,
        "Greywacke": 0.7,
        "Default": 0.6,
    }

    MATERIAL_MARKERS: ClassVar[dict] = {
        "Gneiss": "s",
        "Granite": "o",
        "Greywacke": "^",
        "Default": "D",
    }

    LOAD_LINESTYLES: ClassVar[dict] = {
        "A": "-",
        "B": "--",
        "C": "-.",
        "E": ":",
        "F": (0, (3, 5, 1, 5)),
        "Default": "-",
    }

    def __init__(
        self,
        output_dir,
        colorbar_opts=None,
        markers=None,
        save_extracted_data=False,
    ):
        self.output_dir = Path(output_dir)
        self.save_extracted_data = save_extracted_data
        self.colorbar_opts = colorbar_opts or {
            "u": {"vmin": 0, "vmax": 0.25, "cmap": "Greens"},
            "stress": {"vmin": -20, "vmax": 10, "cmap": "coolwarm"},
            "strain": {"vmin": -0.05, "vmax": 0, "cmap": "RdBu_r"},
            "pressure": {"vmin": 0.1, "vmax": 5, "cmap": "jet"},
        }
        self.load_order = ["A", "B", "C", "D", "E", "F", "NA"]
        self.markers = markers or ["o", "s", "D", "^"]
        self.vtu_file_names: dict[str, dict[str, list[str]]] = {}
        self.material_names: list[str] = []

    def setup_plot_style(self):
        plt.rcParams.update(self.DEFAULT_RC_PARAMS)

    def get_style(self, method=None, material=None, load_case=None):
        base_color = self.METHOD_COLORS.get(method, self.METHOD_COLORS["Unknown"])
        shade = self.MATERIAL_SHADES.get(material, self.MATERIAL_SHADES["Default"])
        color = self._adjust_shade(base_color, shade)
        marker = self.MATERIAL_MARKERS.get(material, self.MATERIAL_MARKERS["Default"])
        linestyle = self.LOAD_LINESTYLES.get(load_case, self.LOAD_LINESTYLES["Default"])
        return color, marker, linestyle

    def get_mesh(self, filename: str, r=0.065, inner=False) -> ot.Mesh:
        mesh = ot.MeshSeries(filename)[-1]
        if inner:
            radii = np.asarray([np.linalg.norm(cell.center) for cell in mesh.cell])
            return ot.Mesh(mesh.extract_cells(radii < r))
        return mesh

    def _plot_mesh(self, filename, save_file=None, r=0.065, inner=False):
        mesh = self.get_mesh(filename, r, inner)
        displacement = ot.variables.displacement.replace(output_unit="mm")

        has_pressure = "pressure" in mesh.point_data
        has_phasefield = "phasefield" in mesh.point_data
        ncols = 4 + int(has_phasefield) if has_pressure else 3 + int(has_phasefield)
        fig, axs = plt.subplots(nrows=1, ncols=ncols, figsize=(12 * ncols, 10))

        axs[0].set_title(r"$\mathbf{u}$ [mm]", fontsize=24)
        mesh.plot_contourf(displacement, fig, axs[0], **self.colorbar_opts["u"])

        axs[1].set_title(r"$\mathrm{tr}(\sigma)$ [MPa]", fontsize=24)
        mesh.plot_contourf(
            ot.variables.stress.trace,
            fig,
            axs[1],
            **self.colorbar_opts["stress"],
        )

        axs[2].set_title(r"$\mathrm{tr}(\varepsilon)$", fontsize=24)
        mesh.plot_contourf(
            ot.variables.strain.trace,
            fig,
            axs[2],
            **self.colorbar_opts["strain"],
        )

        col_idx = 3
        if has_pressure:
            axs[col_idx].set_title(r"$p$ [MPa]", fontsize=24)
            mesh.plot_contourf(
                ot.variables.pressure,
                fig,
                axs[col_idx],
                **self.colorbar_opts["pressure"],
            )
            col_idx += 1

        if has_phasefield:
            var_name = "phasefield" if "phasefield" in mesh.point_data else "d"
            axs[col_idx].set_title(r"$v$ (Phase Field)", fontsize=24)
            mesh.plot_contourf(
                var_name,
                fig,
                axs[col_idx],
                **self.colorbar_opts.get("phasefield", {}),
            )

        domain_limit = 0.098
        for ax in axs:
            ax.set_aspect("equal")
            ax.tick_params(labelsize=16)
            ax.set_xlim(-domain_limit, domain_limit)
            ax.set_ylim(-domain_limit, domain_limit)
        fig.tight_layout()
        if save_file:
            Path(save_file).parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(save_file, dpi=300, bbox_inches="tight")
        plt.show()
        plt.close(fig)

    def _make_filename(self, base_path, suffix):
        base = Path(base_path)
        return base.with_name(f"{base.stem}_{suffix}{base.suffix}")

    def plot_field_variables(self, vtu_input, save_path=None, r=0.065, inner=False):
        def extract_benchmark_and_material(file_path):
            parts = Path(file_path).stem.split("_")
            benchmark = parts[0] if len(parts) > 0 else "Unknown"
            material = parts[-1] if len(parts) > 1 else "Unknown"
            return benchmark, material

        if isinstance(vtu_input, str):
            benchmark, material = extract_benchmark_and_material(vtu_input)
            print(f"Benchmark: {benchmark}, Material: {material}")
            return self._plot_mesh(vtu_input, save_path, r, inner)

        if isinstance(vtu_input, list):
            for i, f in enumerate(vtu_input):
                benchmark, material = extract_benchmark_and_material(f)
                print(f"[{i}] Benchmark: {benchmark}, Material: {material}")
                sp = self._make_filename(save_path, i) if save_path else None
                self._plot_mesh(f, sp, r, inner)
            return None

        if isinstance(vtu_input, dict):
            for key, file_list in vtu_input.items():
                for i, f in enumerate(file_list):
                    benchmark, material = extract_benchmark_and_material(f)
                    tag = f"{key}_{i}" if len(file_list) > 1 else f"{key}"
                    print(
                        f"Load: {key}, File {i}, Benchmark: {benchmark}, Material: {material}"
                    )
                    sp = self._make_filename(save_path, tag) if save_path else None
                    self._plot_mesh(f, sp, r, inner)
            return None

        msg = "vtu_input must be str, list, or dict"
        raise TypeError(msg)

    def inner_mesh(self, filename: str) -> ot.Mesh:
        mesh = ot.MeshSeries(self.output_dir / filename)[-1]
        radii = np.asarray([np.linalg.norm(cell.center) for cell in mesh.cell])
        return ot.Mesh(mesh.extract_cells(radii < 0.065))

    def sorted_angles_eps_trace(self, filename: str) -> tuple[np.ndarray, np.ndarray]:
        edge = self.inner_mesh(filename).extract_feature_edges()
        eps_trace = ot.variables.strain.trace.transform(edge)
        phi = np.arctan2(edge.points[:, 1], edge.points[:, 0])
        phi[phi < 0] += 2 * np.pi
        sort_idx = np.argsort(phi)
        return phi[sort_idx], eps_trace[sort_idx]

    def plot_volumetric_strain_vs_angle(
        self,
        vtu_files,
        model_type="M1",
        ylim_range=None,
        layout="single",
        downsample: int = 1,
        markevery: int | None = None,
        external_data=None,
    ):
        self.setup_plot_style()
        BASE_WIDTH, BASE_HEIGHT = 7.0, 4.0

        if isinstance(vtu_files, list):
            vtu_files = {"Default": vtu_files}
        load_labels = list(vtu_files.keys())

        if layout == "single":
            fig, ax = plt.subplots(1, 1, figsize=(BASE_WIDTH, BASE_HEIGHT))
            axes = dict.fromkeys(load_labels, ax)
        elif layout == "subplots":
            n = len(load_labels)
            fig, axs = plt.subplots(
                n, 1, figsize=(BASE_WIDTH, BASE_HEIGHT * n), sharex=True
            )
            if n == 1:
                axs = [axs]
            axes = {lbl: axs[i] for i, lbl in enumerate(load_labels)}
        else:
            err_msg = f"Invalid layout '{layout}'. Choose 'single' or 'subplots'."
            raise ValueError(err_msg)

        for load_label, file_list in vtu_files.items():
            ax = axes[load_label]
            if layout == "subplots":
                ax.set_title(f"Load case: {load_label}", fontsize=10, pad=10)

            for vtu_file in file_list:
                if not Path(vtu_file).exists():
                    continue
                try:
                    phi, eps_v = self.sorted_angles_eps_trace(vtu_file)
                except Exception:
                    continue

                if self.save_extracted_data:
                    save_path = (
                        Path(self.output_dir or ".")
                        / f"extracted_{Path(vtu_file).stem}_volStrain.npz"
                    )
                    np.savez_compressed(save_path, phi=phi, eps_v=eps_v)

                phi_deg = phi * 180 / math.pi
                eps_pct = np.array(eps_v) * 100
                phi_ds = phi_deg[::downsample]
                eps_ds = eps_pct[::downsample]
                mk = markevery or max(1, len(phi_ds) // 10)

                parts = Path(vtu_file).stem.split("_")
                method = parts[1] if len(parts) > 1 else "Unknown"
                material = parts[3].title() if len(parts) > 3 else "Default"

                base = self.METHOD_COLORS.get(method, self.METHOD_COLORS["Unknown"])
                shade = self.MATERIAL_SHADES.get(
                    material, self.MATERIAL_SHADES["Default"]
                )
                color = self._adjust_shade(base, shade)
                marker = self.MATERIAL_MARKERS.get(
                    material, self.MATERIAL_MARKERS["Default"]
                )
                ls = self.LOAD_LINESTYLES.get(
                    load_label, self.LOAD_LINESTYLES["Default"]
                )
                lbl = f"{material}, {method}, Load {load_label}"

                ax.plot(
                    phi_ds,
                    eps_ds,
                    label=lbl,
                    color=color,
                    linestyle=ls,
                    marker=marker,
                    markersize=6,
                    markevery=mk,
                    markerfacecolor="none",
                    markeredgewidth=1.2,
                )

            if external_data:
                for method_name, load_dict in external_data.items():
                    if load_label not in load_dict:
                        continue
                    for mat_name, vals in load_dict[load_label].items():
                        try:
                            phi_ext = np.array(vals["angle"])
                            eps_ext = np.array(vals["vol_strain"]) * 100
                        except Exception:
                            continue

                        phi_ds_ext = phi_ext[::downsample]
                        eps_ds_ext = eps_ext[::downsample]
                        mk_ext = markevery or max(1, len(phi_ds_ext) // 10)

                        base = self.METHOD_COLORS.get(
                            method_name, self.METHOD_COLORS["Unknown"]
                        )
                        shade = self.MATERIAL_SHADES.get(
                            mat_name.title(), self.MATERIAL_SHADES["Default"]
                        )
                        color = self._adjust_shade(base, shade)
                        marker = self.MATERIAL_MARKERS.get(
                            mat_name.title(), self.MATERIAL_MARKERS["Default"]
                        )
                        ls = self.LOAD_LINESTYLES.get(
                            load_label, self.LOAD_LINESTYLES["Default"]
                        )
                        lbl = f"{mat_name.title()}, {method_name}, Load {load_label}"

                        ax.plot(
                            phi_ds_ext,
                            eps_ds_ext,
                            label=lbl,
                            color=color,
                            linestyle=ls,
                            marker=marker,
                            markersize=6,
                            markevery=mk_ext,
                            markerfacecolor="none",
                            markeredgewidth=1.2,
                        )

            ax.set_xlim(0, 360)
            ax.set_xticks(np.arange(0, 361, 45))
            ax.grid(True, which="major", linestyle=":", linewidth=0.5, alpha=0.7)
            ax.set_xlabel(r"Angle $\theta$ [$^\circ$]", labelpad=5)
            ax.set_ylabel(
                r"Volumetric strain $\varepsilon_{\mathrm{vol}}$ [%]",
                labelpad=5,
            )
            ax.set_ylim(ylim_range or [-10, 0])
            ax.minorticks_on()
            ax.tick_params(which="minor", length=3, width=0.6)
            ax.tick_params(which="major", length=5, width=0.8)
            ax.tick_params(
                axis="both", which="both", direction="in", top=True, right=True
            )

        Path(self.output_dir or ".").mkdir(parents=True, exist_ok=True)
        save_opts = {"dpi": 600, "bbox_inches": "tight", "pad_inches": 0.05}

        if layout == "single":
            handles, labels = ax.get_legend_handles_labels()
            ax.legend(
                handles,
                labels,
                loc="lower right",
                bbox_to_anchor=(0.98, 0.2),
                frameon=True,
                framealpha=0.95,
                edgecolor="0.3",
                borderpad=0.6,
            )
            fig.tight_layout()
            fig.savefig(
                Path(self.output_dir, f"volume_strain_plot_{model_type}_all.pdf"),
                **save_opts,
            )
        else:
            fig.subplots_adjust(hspace=0.3)
            for subax in fig.axes:
                handles, labels = subax.get_legend_handles_labels()
                subax.legend(
                    handles,
                    labels,
                    loc="lower right",
                    bbox_to_anchor=(0.98, 0.098),
                    frameon=True,
                    framealpha=0.95,
                    edgecolor="0.3",
                    borderpad=0.6,
                )
            fig.tight_layout()
            fig.savefig(
                Path(
                    self.output_dir,
                    f"volume_strain_plot_{model_type}_subplots.pdf",
                ),
                **save_opts,
            )

        plt.show()
        plt.close()

    def _adjust_shade(self, hex_color, factor):
        try:
            hex_color = hex_color.strip("#")
            if len(hex_color) != 6:
                msg = f"Invalid hex color: #{hex_color}"
                raise ValueError(msg)

            rgb = tuple(int(hex_color[i : i + 2], 16) for i in (0, 2, 4))
            h, lum, s = colorsys.rgb_to_hls(*(x / 255.0 for x in rgb))
            lum = max(0.2, min(0.9, factor * (lum or 0.5)))
            rgb = colorsys.hls_to_rgb(h, lum, s)
            return "#{:02x}{:02x}{:02x}".format(*tuple(int(x * 255) for x in rgb))
        except Exception as e:
            print(f"Error adjusting shade for color {hex_color}: {e}")
            return "#000000"

    def get_sub_mesh_by_material(self, mesh: ot.Mesh, material_id: int) -> ot.Mesh:
        if "MaterialIDs" not in mesh.cell_data:
            msg = f"Mesh missing 'MaterialIDs' in cell_data; cannot extract material {material_id}."
            raise KeyError(msg)

        mat_ids = mesh.cell_data["MaterialIDs"]
        idx = np.where(mat_ids == material_id)[0]
        return mesh.extract_cells(idx)

    def _avg_fracture_aperture(self, vtu_file: str, material_id: int) -> float:
        mesh = ot.MeshSeries(self.output_dir / vtu_file)[-1]
        sub = self.get_sub_mesh_by_material(mesh, material_id)
        if "fracture_aperture" not in sub.point_data:
            return np.nan
        aperture_vals = sub.point_data["fracture_aperture"]
        valid = aperture_vals[aperture_vals >= 0]
        return float(np.mean(valid)) if valid.size > 0 else np.nan

    def _avg_fracture_normal_stress(self, vtu_file: str, material_id: int) -> float:
        mesh = ot.MeshSeries(self.output_dir / vtu_file)[-1]
        sub = self.get_sub_mesh_by_material(mesh, material_id)
        if "fracture_stress" not in sub.point_data:
            return np.nan
        stress_all = sub.point_data["fracture_stress"]
        normal_vals = stress_all[:, 1]
        return float(np.mean(normal_vals)) if normal_vals.size > 0 else np.nan

    def extract_lie_aperture_from_list(
        self, file_list: list[str], material_names: list[str], fracture_mat_id=3
    ) -> dict[str, tuple[np.ndarray, np.ndarray]]:
        by_material: dict[str, tuple[np.ndarray, np.ndarray]] = {}
        for vtu_file, mat_name in zip(file_list, material_names, strict=True):
            mesh = ot.MeshSeries(self.output_dir / vtu_file)[-1]
            sub = self.get_sub_mesh_by_material(mesh, fracture_mat_id)
            coords = sub.points
            vals = sub.point_data.get("fracture_aperture", np.array([]))
            by_material[mat_name] = (coords, vals)
        return by_material

    def _gather_lie_avg_data(self) -> tuple[list[str], np.ndarray, np.ndarray]:
        lie_dict = self.vtu_file_names.get("LIE", {})
        load_cases = list(lie_dict.keys())
        n_loads = len(load_cases)
        n_mat = len(self.material_names)
        avg_stress_n_lie = np.full((n_loads, n_mat), np.nan, dtype=float)
        avg_aperture_lie = np.full((n_loads, n_mat), np.nan, dtype=float)
        for i, load_case in enumerate(load_cases):
            file_list = lie_dict[load_case]
            for j, _ in enumerate(self.material_names):
                vtu_file = file_list[j]
                avg_stress = self._avg_fracture_normal_stress(vtu_file, material_id=3)
                avg_aperture = self._avg_fracture_aperture(vtu_file, material_id=3)
                avg_stress_n_lie[i, j] = avg_stress
                avg_aperture_lie[i, j] = avg_aperture
        return load_cases, avg_stress_n_lie, avg_aperture_lie

    def extract_vpf_width_from_list(
        self, file_list: list[str], material_names: list[str]
    ) -> dict[str, tuple[np.ndarray, np.ndarray]]:
        width_by_material: dict[str, tuple[np.ndarray, np.ndarray]] = {}
        for vtu_file, material_name in zip(file_list, material_names, strict=True):
            mesh = ot.MeshSeries(self.output_dir / vtu_file)[-1]
            cell_data = mesh.cell_data

            if "damage" not in cell_data or "width" not in cell_data:
                width_by_material[material_name] = None
                continue

            damage = cell_data["damage"]
            width = cell_data["width"]
            damaged_idx = np.where(damage == 0)[0]

            if damaged_idx.size == 0:
                width_by_material[material_name] = None
                continue

            w_vals = np.abs(width[damaged_idx])
            coords = mesh.cell_centers().points[damaged_idx]
            width_by_material[material_name] = (coords, w_vals)

        return width_by_material

    def extract_vpf_stress_from_list(
        self,
        file_list: list[str],
        material_names: list[str],
        stress_type: str = "normal",
    ) -> dict[str, tuple[np.ndarray, np.ndarray]]:
        stress_by_material: dict[str, tuple[np.ndarray, np.ndarray]] = {}
        for vtu_file, material_name in zip(file_list, material_names, strict=True):
            mesh = ot.MeshSeries(self.output_dir / vtu_file)[-1]
            cell_data = mesh.cell_data
            point_data = mesh.point_data

            if "damage" not in cell_data or "sigma" not in point_data:
                stress_by_material[material_name] = None
                continue

            damage = cell_data["damage"]
            sigma = point_data["sigma"]
            damaged_idx = np.where(damage == 0)[0]

            if damaged_idx.size == 0:
                stress_by_material[material_name] = None
                continue

            cells_tri = mesh.extract_cells(mesh.celltypes == 5)
            cells_quad = mesh.extract_cells(mesh.celltypes == 9)
            tri_conn = (
                cells_tri.cells.reshape(-1, 4)[:, 1:]
                if cells_tri.n_cells > 0
                else np.empty((0, 3), dtype=int)
            )
            quad_conn = (
                cells_quad.cells.reshape(-1, 5)[:, 1:]
                if cells_quad.n_cells > 0
                else np.empty((0, 4), dtype=int)
            )
            coords_list, stress_list = [], []
            centers = mesh.cell_centers().points

            for cidx in damaged_idx:
                if cidx >= tri_conn.shape[0] + quad_conn.shape[0]:
                    continue

                if cidx < tri_conn.shape[0]:
                    node_ids = tri_conn[cidx]
                else:
                    node_ids = quad_conn[cidx - tri_conn.shape[0]]

                coords_list.append(centers[cidx])
                node_sigma = sigma[node_ids]

                if stress_type == "von_mises":
                    sigma_xx = node_sigma[:, 0]
                    sigma_yy = node_sigma[:, 1]
                    tau_xy = node_sigma[:, 3]
                    vm = np.sqrt(
                        sigma_xx**2
                        + sigma_yy**2
                        - sigma_xx * sigma_yy
                        + 3 * tau_xy**2
                    )
                    stress_list.append(np.mean(vm))
                else:
                    sigma_yy = node_sigma[:, 1]
                    stress_list.append(np.mean(sigma_yy))

            if not coords_list:
                stress_by_material[material_name] = None
                continue

            stress_by_material[material_name] = (
                np.vstack(coords_list),
                np.array(stress_list),
            )
        return stress_by_material

    def get_vpf_permeability_data(self) -> tuple[np.ndarray, np.ndarray]:
        vpf_dict = self.vtu_file_names.get("VPF", {})
        mean_w2_list, mean_w_list = [], []
        for _, file_list in vpf_dict.items():
            w2_row, w_row = [], []
            for vtu_file in file_list:
                mesh = ot.MeshSeries(self.output_dir / vtu_file)[-1]
                cell_data = mesh.cell_data

                if "damage" not in cell_data or "width" not in cell_data:
                    w2_row.append(np.nan)
                    w_row.append(np.nan)
                    continue

                damage = cell_data["damage"]
                width = cell_data["width"]
                damaged_idx = np.where(damage == 0)[0]

                if damaged_idx.size == 0:
                    w2_row.append(np.nan)
                    w_row.append(np.nan)
                    continue

                w_vals = width[damaged_idx]
                w2_row.append(float(np.mean(w_vals**2)))
                w_row.append(float(np.mean(np.abs(w_vals))))

            mean_w2_list.append(w2_row)
            mean_w_list.append(w_row)
        return np.array(mean_w2_list, dtype=float), np.array(mean_w_list, dtype=float)

    def _gather_vpf_avg_data(self) -> tuple[list[str], np.ndarray, np.ndarray]:
        vpf_dict = self.vtu_file_names.get("VPF", {})
        load_cases = list(vpf_dict.keys())
        n_loads = len(load_cases)
        n_mat = len(self.material_names)
        avg_stress_n_vpf = np.full((n_loads, n_mat), np.nan, dtype=float)
        avg_width_vpf = np.full((n_loads, n_mat), np.nan, dtype=float)
        for i, load_case in enumerate(load_cases):
            file_list = vpf_dict[load_case]
            stress_dict = self.extract_vpf_stress_from_list(
                file_list, self.material_names, stress_type="normal"
            )
            width_dict = self.extract_vpf_width_from_list(
                file_list, self.material_names
            )
            for j, mat_name in enumerate(self.material_names):
                s_entry = stress_dict.get(mat_name)
                if s_entry is not None:
                    _, s_vals = s_entry
                    if s_vals is not None and s_vals.size > 0:
                        avg_stress_n_vpf[i, j] = np.nanmean(s_vals)
                w_entry = width_dict.get(mat_name)
                if w_entry is not None:
                    _, w_vals = w_entry
                    if w_vals is not None and w_vals.size > 0:
                        avg_width_vpf[i, j] = np.nanmean(w_vals)
        return load_cases, avg_stress_n_vpf, avg_width_vpf

    def compute_f_stress_from_pee(
        self, pee_load_values: dict[str, list[float]], step_deg: float = 22.5
    ) -> dict[str, float]:
        step_rad = math.radians(step_deg)
        result: dict[str, float] = {}
        for case, values in pee_load_values.items():
            if not values:
                result[case] = float("nan")
                continue
            s2 = min(values)
            s3 = max(values)
            idx3 = values.index(s3)
            theta = idx3 * step_rad
            f_val = abs(0.5 * (s2 + s3) + 0.5 * (s2 - s3) * math.cos(2 * theta))
            result[case] = f_val
        return result

    def plot_avg_width_vs_stress(
        self,
        pee_load_values: dict[str, list[float]],
        metric: str = "width",
        methods_to_include: list[str] | None = None,
        ylim_range: tuple[float, float] | None = None,
        external_data: dict[str, dict[str, dict[str, list[float]]]] | None = None,
        benchmark_tag: str | None = None,
    ) -> np.ndarray:
        if methods_to_include is None:
            methods_to_include = ["VPF", "LIE"]
        self.setup_plot_style()

        comp = self.compute_f_stress_from_pee(pee_load_values)
        f_arr = np.array([comp[lc] for lc in pee_load_values])
        self.f_stress_n_estimated = f_arr
        load_cases = list(pee_load_values.keys())

        lc_vpf, stress_vpf, width_vpf = ([], np.empty((0, 0)), np.empty((0, 0)))
        lc_lie, stress_lie, aperture_lie = ([], np.empty((0, 0)), np.empty((0, 0)))
        if "VPF" in methods_to_include:
            lc_vpf, stress_vpf, width_vpf = self._gather_vpf_avg_data()
        if "LIE" in methods_to_include:
            lc_lie, stress_lie, aperture_lie = self._gather_lie_avg_data()

        base = lc_vpf or lc_lie
        idx_map = {case: i for i, case in enumerate(base)}
        indices = [idx_map[lc] for lc in load_cases]
        if stress_vpf.size:
            stress_vpf = stress_vpf[indices, :]
            width_vpf = width_vpf[indices, :]
        if stress_lie.size:
            stress_lie = stress_lie[indices, :]
            aperture_lie = aperture_lie[indices, :]

        method_data: dict[str, dict[str, np.ndarray]] = {}
        if stress_vpf.size:
            arr = {
                "width": width_vpf,
                "stress": stress_vpf,
                "permeability": width_vpf**2 / 12.0,
            }[metric]
            method_data["VPF"] = {
                mat: arr[:, j] for j, mat in enumerate(self.material_names)
            }
        if stress_lie.size:
            arr = {
                "width": aperture_lie,
                "stress": stress_lie,
                "permeability": aperture_lie**2 / 12.0,
            }[metric]
            method_data["LIE"] = {
                mat: arr[:, j] for j, mat in enumerate(self.material_names)
            }

        fig, ax = plt.subplots(figsize=(7, 4.5))
        for method, data in method_data.items():
            base_color = self.METHOD_COLORS.get(method)
            for mat, y in data.items():
                shade = self.MATERIAL_SHADES.get(mat, self.MATERIAL_SHADES["Default"])
                color = self._adjust_shade(base_color, shade)
                marker = self.MATERIAL_MARKERS.get(mat)
                sample_load = load_cases[list(data.keys()).index(mat) % len(load_cases)]
                ls = self.LOAD_LINESTYLES.get(sample_load)
                ax.plot(
                    f_arr,
                    y,
                    color=color,
                    linestyle=ls,
                    marker=marker,
                    markersize=6,
                    markevery=1,
                    markerfacecolor="none",
                    markeredgewidth=1.2,
                    label=f"{mat}, {method}",
                )

        first_method = next(iter(method_data), None)
        if first_method:
            first_mat = next(iter(method_data[first_method]))
            y0 = method_data[first_method][first_mat]
            for idx, lc in enumerate(load_cases):
                ax.annotate(
                    lc,
                    (f_arr[idx], y0[idx]),
                    textcoords="offset points",
                    xytext=(0, 8),
                    ha="center",
                    fontsize=8,
                )

        if external_data:
            for method_name, load_dict in external_data.items():
                base_color = self.METHOD_COLORS.get(method_name)
                for load_case, mat_profiles in load_dict.items():
                    if load_case != "all" and load_case not in load_cases:
                        continue
                    ls = self.LOAD_LINESTYLES.get(load_case, "-")
                    for mat, vals in mat_profiles.items():
                        s_ext = np.array(vals.get("f_stress", []))
                        y_ext = np.array(
                            vals.get(metric, [])
                            if metric != "permeability"
                            else vals.get("permeability", [])
                        )
                        if s_ext.size and y_ext.size:
                            shade = self.MATERIAL_SHADES.get(
                                mat, self.MATERIAL_SHADES["Default"]
                            )
                            color_ext = self._adjust_shade(base_color, shade)
                            marker_ext = self.MATERIAL_MARKERS.get(mat)
                            ax.plot(
                                s_ext,
                                y_ext,
                                color=color_ext,
                                linestyle=ls,
                                marker=marker_ext,
                                markersize=6,
                                markevery=1,
                                markerfacecolor="none",
                                markeredgewidth=1.2,
                                label=f"{mat}, {method_name}",
                            )

        ax.set_xlabel("Fracture normal stress [Pa]", labelpad=5)
        ylabel_map = {
            "width": "Average width [m]",
            "stress": "Average stress [MPa]",
            "permeability": "Permeability [m²]",
        }
        ax.set_ylabel(ylabel_map[metric], labelpad=5)
        if ylim_range:
            ax.set_ylim(ylim_range)
        ax.grid(True, which="major", linestyle=":", linewidth=0.5, alpha=0.7)
        ax.minorticks_on()
        ax.tick_params(which="minor", length=3, width=0.6)
        ax.tick_params(which="major", length=5, width=0.8)
        ax.tick_params(axis="both", which="both", direction="in", top=True, right=True)

        handles, labels = ax.get_legend_handles_labels()
        ax.legend(
            handles,
            labels,
            loc="upper right",
            frameon=True,
            framealpha=0.95,
            edgecolor="0.3",
            borderpad=0.6,
        )

        if self.save_extracted_data:
            tag = benchmark_tag or "Unknown"
            for method, data in method_data.items():
                for mat, y in data.items():
                    outp = (
                        self.output_dir
                        / f"extracted_{tag}_{method}_{mat}_avg_{metric}.npz"
                    )
                    np.savez_compressed(outp, f_stress=f_arr, **{metric: y})

        if benchmark_tag:
            out_pdf = (
                self.output_dir / f"avg_width_vs_stress_{benchmark_tag}_{metric}.pdf"
            )
            out_pdf.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(out_pdf, dpi=600, bbox_inches="tight", pad_inches=0.05)

        plt.tight_layout()
        plt.show()
        plt.close(fig)

        return f_arr

    def plot_fracture_aperture_profiles(
        self,
        widthProfile: dict[str, dict[str, tuple[np.ndarray, np.ndarray]]],
        benchmark_tag: str,
        downsample: int = 1,
        markevery: int | None = None,
        ylim: tuple[float, float] | None = None,
        method_label: str = "FEM",
        external_data: dict[str, dict[str, dict[str, dict[str, list[float]]]]]
        | None = None,
    ) -> None:
        self.setup_plot_style()
        for load_case, mat_prof in widthProfile.items():
            fig, ax = plt.subplots(figsize=(8, 5))
            ax.set_title(f"Load case: {load_case}", fontsize=10, pad=10)
            for mat, (coords_arr, vals_arr) in mat_prof.items():
                if (
                    coords_arr is None
                    or vals_arr is None
                    or coords_arr.size == 0
                    or vals_arr.size == 0
                ):
                    continue
                idx = np.argsort(coords_arr[:, 0])
                coords, vals = coords_arr[idx], vals_arr[idx]

                if self.save_extracted_data:
                    outp = (
                        self.output_dir
                        / f"extracted_{benchmark_tag}_{method_label}_{load_case}_{mat}_profile.npz"
                    )
                    np.savez_compressed(outp, coords=coords, aperture=vals)

                mk = markevery or max(1, len(coords[::downsample, 0]) // 10)
                base = self.METHOD_COLORS.get(
                    method_label, self.METHOD_COLORS["Default"]
                )
                shade = self.MATERIAL_SHADES.get(mat, self.MATERIAL_SHADES["Default"])
                color = self._adjust_shade(base, shade)
                marker = self.MATERIAL_MARKERS.get(
                    mat, self.MATERIAL_MARKERS["Default"]
                )
                ls = self.LOAD_LINESTYLES.get(load_case, "-")
                lbl = f"{mat}, {method_label}"

                ax.plot(
                    coords[::downsample, 0],
                    vals[::downsample],
                    color=color,
                    linestyle=ls,
                    marker=marker,
                    markersize=6,
                    markevery=mk,
                    markerfacecolor="none",
                    markeredgewidth=1.2,
                    label=lbl,
                )

            if external_data:
                for ext_method, method_dict in external_data.items():
                    if load_case not in method_dict:
                        continue
                    for mat, rec in method_dict[load_case].items():
                        coords_ext = np.array(rec.get("coords", []))
                        aper_ext = np.array(rec.get("aperture", []))
                        if coords_ext.size == 0 or aper_ext.size == 0:
                            continue
                        idx_ext = np.argsort(coords_ext[:, 0])
                        coords_ext, aper_ext = (
                            coords_ext[idx_ext],
                            aper_ext[idx_ext],
                        )

                        if self.save_extracted_data:
                            outp = (
                                self.output_dir
                                / f"extracted_{benchmark_tag}_{ext_method}_{load_case}_{mat}_profile.npz"
                            )
                            np.savez_compressed(
                                outp, coords=coords_ext, aperture=aper_ext
                            )

                        mk_ext = markevery or max(
                            1, len(coords_ext[::downsample, 0]) // 10
                        )

                        base = self.METHOD_COLORS.get(
                            ext_method, self.METHOD_COLORS["Default"]
                        )
                        shade = self.MATERIAL_SHADES.get(
                            mat, self.MATERIAL_SHADES["Default"]
                        )
                        color = self._adjust_shade(base, shade)
                        marker = self.MATERIAL_MARKERS.get(
                            mat, self.MATERIAL_MARKERS["Default"]
                        )
                        ls = self.LOAD_LINESTYLES.get(load_case, "-")
                        lbl_ext = f"{mat}, {ext_method}"

                        ax.plot(
                            coords_ext[::downsample, 0],
                            aper_ext[::downsample],
                            color=color,
                            linestyle=ls,
                            marker=marker,
                            markersize=6,
                            markevery=mk_ext,
                            markerfacecolor="none",
                            markeredgewidth=1.2,
                            label=lbl_ext,
                        )

            ax.set_xlabel("x [m]", labelpad=8)
            ax.set_ylabel("Aperture [m]", labelpad=8)
            if ylim is None:
                ymax = max(
                    np.nanmax(vals) * 1.1
                    for _, vals in mat_prof.values()
                    if vals is not None and vals.size > 0
                )
                ylim = (0, max(1e-6, ymax))
            ax.set_ylim(ylim)
            if ylim[1] < 1e-3:
                ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
            ax.grid(True, which="both", linestyle=":", alpha=0.5)
            ax.legend(
                loc="upper right",
                frameon=True,
                framealpha=0.95,
                edgecolor="0.8",
            )

            out_file = (
                self.output_dir / f"aperture_profiles_{load_case}_{benchmark_tag}.pdf"
            )
            fig.savefig(out_file, dpi=600, bbox_inches="tight", pad_inches=0.05)
            plt.tight_layout()
            plt.show()
            plt.close()

    @staticmethod
    def load_external_data(
        ext_dir: str | Path, benchmark_tag: str | None = None
    ) -> dict[str, dict]:
        """
        Scans ext_dir for extracted_*.npz and returns dicts:
          - "strain": volumetric-strain data
          - "average": avg-data (width/stress/permeability)
          - "widthProfile": profile data
        """
        ext = Path(ext_dir)
        raw: dict[str, dict] = {}

        for p in ext.glob("extracted_*.npz"):
            parts = p.stem.split("_")
            if parts[0] != "extracted":
                continue

            # volumetric-strain
            if len(parts) == 6 and parts[-1] == "volStrain":
                _, bench, method, load, mat, _ = parts
                suffix = "volStrain"

            # average
            elif (len(parts) == 7 and parts[5] == "avg") or (
                len(parts) == 6 and parts[4] == "avg"
            ):
                if len(parts) == 7:
                    _, bench, method, load, mat, _, metric = parts
                else:
                    _, bench, method, mat, _, metric = parts
                    load = "all"
                suffix = f"avg_{metric}"

            # widthProfile
            elif len(parts) >= 6 and parts[-1] == "profile":
                _, bench, method, load, mat, *_ = parts
                suffix = "profile"
            else:
                continue

            if benchmark_tag is not None and bench != benchmark_tag:
                continue

            arr = np.load(p)
            data = {k: arr[k].tolist() for k in arr.files}
            raw.setdefault(method, {}).setdefault(load, {}).setdefault(mat, {})[
                suffix
            ] = data

        strain: dict = {}
        average: dict = {}
        widthProfile: dict = {}

        for method, loads in raw.items():
            for load, mats in loads.items():
                for mat, sufmap in mats.items():
                    # volumetric strain
                    if "volStrain" in sufmap:
                        sd = sufmap["volStrain"]
                        strain.setdefault(method, {}).setdefault(load, {})[mat] = {
                            "angle": [phi * 180 / np.pi for phi in sd["phi"]],
                            "vol_strain": sd["eps_v"],
                        }
                    # average f_stress + metric
                    for key, sd in sufmap.items():
                        if key.startswith("avg_"):
                            metric = key.split("_", 1)[1]
                            average.setdefault(method, {}).setdefault(
                                load, {}
                            ).setdefault(mat, {})["f_stress"] = sd["f_stress"]
                            average[method][load][mat][metric] = sd[metric]
                    # profile
                    if "profile" in sufmap:
                        pd = sufmap["profile"]
                        widthProfile.setdefault(method, {}).setdefault(load, {})[
                            mat
                        ] = {
                            "coords": pd["coords"],
                            "aperture": pd["aperture"],
                        }

        def _sort(d):
            if not isinstance(d, dict):
                return d
            return {k: _sort(d[k]) for k in sorted(d)}

        return {
            "strain": _sort(strain),
            "average": _sort(average),
            "widthProfile": _sort(widthProfile),
        }

    def plot_observation_points_vs_time(
        self,
        vtu_files,
        obs_points,
        var,
        methods_to_include=None,
        materials_to_include=None,
        markevery=None,
        time_min=None,
        time_max=None,
        **plot_kwargs,
    ):
        self.setup_plot_style()
        fig, ax = plt.subplots(1, 1, figsize=(7, 4.5))

        for load_label, file_list in vtu_files.items():
            for pvd in file_list:
                parts = Path(pvd).stem.split("_")
                benchmark = parts[0] if parts else "Unknown"
                method = parts[1] if len(parts) > 1 else "Unknown"
                material = parts[-1] if parts else "Default"

                if methods_to_include and method not in methods_to_include:
                    continue
                if materials_to_include and material not in materials_to_include:
                    continue

                ms = ot.MeshSeries(self.output_dir / pvd).scale(time=("s", "s"))
                for obs_label, pts in obs_points.items():
                    probe = ot.MeshSeries.extract_probe(ms, pts)
                    times = probe.timevalues
                    values = var.transform(probe)

                    mask = np.ones_like(times, dtype=bool)
                    if time_min is not None:
                        mask &= times >= time_min
                    if time_max is not None:
                        mask &= times <= time_max
                    times, values = times[mask], values[mask, :]

                    if self.save_extracted_data:
                        # Save one file per observation series
                        for pt_idx in range(values.shape[1]):
                            key_file = (
                                f"extracted_{benchmark}_{method}_{load_label}_"
                                f"{material}_{obs_label}_pt{pt_idx}.npz"
                            )
                            out_path = self.output_dir / key_file
                            np.savez_compressed(
                                out_path,
                                times=times,
                                values=values[:, pt_idx],
                                points=np.array(pts[pt_idx]),
                            )

                    base_col = self.METHOD_COLORS.get(method)
                    shade = self.MATERIAL_SHADES.get(material)
                    color = self._adjust_shade(base_col, shade)
                    marker = self.MATERIAL_MARKERS.get(material)
                    ls = self.LOAD_LINESTYLES.get(load_label)

                    for pt_idx in range(values.shape[1]):
                        coord_str = np.array2string(
                            np.asarray(pts[pt_idx]),
                            precision=3,
                            separator=", ",
                            suppress_small=True,
                            floatmode="fixed",
                        )
                        lbl = f"{method}, {material}, {load_label}, {obs_label} {coord_str}"
                        ax.plot(
                            times,
                            values[:, pt_idx],
                            label=lbl,
                            color=color,
                            linestyle=ls,
                            marker=marker,
                            markersize=6,
                            markevery=markevery,
                            markerfacecolor="none",
                            markeredgewidth=1.2,
                            **plot_kwargs,
                        )

        ax.set_xlabel("Time [s]", labelpad=5)
        ax.set_ylabel("Pressure [MPa]", labelpad=5)
        ax.grid(True, which="major", linestyle=":", linewidth=0.5, alpha=0.7)
        ax.minorticks_on()
        ax.tick_params(which="minor", length=3, width=0.6)
        ax.tick_params(
            which="major", length=5, width=0.8, direction="in", top=True, right=True
        )

        if time_min is not None or time_max is not None:
            ax.set_xlim(left=time_min, right=time_max)

        handles, labels = ax.get_legend_handles_labels()
        ax.legend(
            handles,
            labels,
            loc="upper right",
            frameon=True,
            framealpha=0.95,
            edgecolor="0.3",
            borderpad=0.6,
        )
        fig.tight_layout()
        out_pdf = self.output_dir / "pressure_obs_vs_time_all.pdf"
        fig.savefig(out_pdf, bbox_inches="tight", pad_inches=0.05)
        plt.show()
        plt.close(fig)
