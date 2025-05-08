import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import ogstools as ot
from matplotlib import cm as colormaps


class Plotter:
    def __init__(
        self,
        output_dir,
        colorbar_opts=None,
        markers=None,
        material_cmaps=None,
        save_extracted_data=False,
    ):
        self.output_dir = output_dir
        self.save_extracted_data = save_extracted_data
        self.colorbar_opts = colorbar_opts or {
            "u": {"vmin": 0, "vmax": 0.25, "cmap": "Greens"},
            "stress": {"vmin": -20, "vmax": 10, "cmap": "coolwarm"},
            "strain": {"vmin": -0.05, "vmax": 0, "cmap": "RdBu_r"},
        }
        self.load_order = ["A", "B", "C", "D", "E", "F", "NA"]

        self.markers = markers or ["o", "s", "D", "^"]
        self.material_cmaps = material_cmaps or {
            "Gneiss": colormaps.get_cmap("Blues").resampled(len(self.load_order)),
            "Greywacke": colormaps.get_cmap("Oranges").resampled(len(self.load_order)),
        }

    def get_mesh(self, filename: str, r=0.065, inner=False) -> ot.Mesh:
        mesh = ot.MeshSeries(Path(self.output_dir, filename))[-1]
        if inner:
            radii = np.asarray([np.linalg.norm(cell.center) for cell in mesh.cell])
            return ot.Mesh(mesh.extract_cells(radii < r))
        return mesh

    def _plot_mesh(self, filename, save_file=None, r=0.065, inner=False):
        mesh = self.get_mesh(filename, r, inner)
        displacement = ot.variables.displacement.replace(output_unit="mm")

        fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(36, 10))

        axs[0].set_title(r"$\mathbf{u}$ [mm]", fontsize=24)
        mesh.plot_contourf(displacement, fig, axs[0], **self.colorbar_opts["u"])

        axs[1].set_title(r"$\mathrm{tr}(\sigma)$ [MPa]", fontsize=24)
        mesh.plot_contourf(
            ot.variables.stress.trace, fig, axs[1], **self.colorbar_opts["stress"]
        )

        axs[2].set_title(r"$\mathrm{tr}(\varepsilon)$", fontsize=24)
        mesh.plot_contourf(
            ot.variables.strain.trace, fig, axs[2], **self.colorbar_opts["strain"]
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
            plt.close(fig)
        else:
            plt.show()

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
            self._plot_mesh(vtu_input, save_path, r, inner)

        elif isinstance(vtu_input, list):
            for i, f in enumerate(vtu_input):
                benchmark, material = extract_benchmark_and_material(f)
                print(f"[{i}] Benchmark: {benchmark}, Material: {material}")
                sp = self._make_filename(save_path, i) if save_path else None
                self._plot_mesh(f, sp, r, inner)

        elif isinstance(vtu_input, dict):
            for key, file_list in vtu_input.items():
                for i, f in enumerate(file_list):
                    benchmark, material = extract_benchmark_and_material(f)
                    tag = f"{key}_{i}" if len(file_list) > 1 else f"{key}"
                    print(
                        f"Load: {key}, File {i}, Benchmark: {benchmark}, Material: {material}"
                    )
                    sp = self._make_filename(save_path, tag) if save_path else None
                    self._plot_mesh(f, sp, r, inner)
        else:
            msg = "vtu_input must be str, list, or dict"
            raise TypeError(msg)

    def inner_mesh(self, filename: str) -> ot.Mesh:
        mesh = ot.MeshSeries(Path(self.output_dir, filename))[-1]
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
        external_data=None,
    ):
        if isinstance(vtu_files, list):
            vtu_files = {"Default": vtu_files}

        load_labels = list(vtu_files.keys())
        all_markers = ["v", "P", "*", "X", "h", "H", "+", "x", "|", "_", ">", "<"]
        external_markers = [m for m in all_markers if m not in self.markers]

        if layout == "single":
            fig, ax = plt.subplots(1, 1, figsize=(12, 4))
            axes = dict.fromkeys(load_labels, ax)
        elif layout == "subplots":
            n = len(load_labels)
            fig, axs = plt.subplots(n, 1, figsize=(12, 4 * n), sharex=True)
            if n == 1:
                axs = [axs]
            axes = {load: axs[i] for i, load in enumerate(load_labels)}

        else:
            msg = f"Invalid layout '{layout}'. Choose from 'single' or 'subplots'."
            raise ValueError(msg)

        for load_label, file_list in vtu_files.items():
            ax = axes[load_label]

            if layout == "subplots":
                ax.set_title(f"Load case: {load_label}", fontsize=14, pad=20)
            # ---- Plot from VTU files ----
            for file_idx, vtu_file in enumerate(file_list):
                if not Path(vtu_file).exists():
                    continue

                try:
                    phi, eps_v = self.sorted_angles_eps_trace(vtu_file)
                except Exception:
                    continue

                if getattr(self, "save_extracted_data", False):
                    save_path = (
                        Path(self.output_dir or ".")
                        / f"extracted_{Path(vtu_file).stem}.npz"
                    )
                    np.savez(save_path, phi=phi, eps_v=eps_v)

                filename = Path(vtu_file).name.replace(".pvd", "")
                parts = filename.split("_")

                method = parts[1] if len(parts) > 1 else "Unknown"
                material = parts[3] if len(parts) > 3 else parts[-1]
                load_letter = load_label

                label = f"{material}, UFZ-OGS {method}, Load {load_letter}"
                material = material.title()

                cmap = self.material_cmaps.get(
                    material, colormaps.get_cmap("gray").resampled(len(self.load_order))
                )
                try:
                    color_idx = self.load_order.index(load_letter)
                except ValueError:
                    color_idx = 0
                color = cmap(color_idx)

                linestyle = "--" if material == "Gneiss" else "-."
                marker = self.markers[file_idx % len(self.markers)]

                ax.plot(
                    phi[::5] * 180 / math.pi,
                    np.array(eps_v[::5]),
                    label=label,
                    color=color,
                    linestyle=linestyle,
                    linewidth=1.8,
                    marker=marker,
                    markersize=6,
                    markerfacecolor="none",
                )

            # ---- Plot from External Data ----
            if external_data:
                for method_name, load_dict in external_data.items():
                    if load_label not in load_dict:
                        continue
                    for material_name, values in load_dict[load_label].items():
                        material = material_name.title()
                        phi = np.array(values["angle"])
                        eps_v = np.array(values["vol_strain"])
                        label = f"{material}, {method_name}, Load {load_label}"

                        cmap = self.material_cmaps.get(
                            material,
                            colormaps.get_cmap("gray").resampled(len(self.load_order)),
                        )
                        try:
                            color_idx = self.load_order.index(load_label)
                        except ValueError:
                            color_idx = 0
                        color = cmap(color_idx)

                        if method_name == "KAERI-TOUGH-3DEC":
                            phi_ds, eps_v_ds = phi, eps_v
                        else:
                            phi_ds = phi[::5]
                            eps_v_ds = eps_v[::5]

                        method_idx = list(external_data.keys()).index(method_name)
                        marker = external_markers[method_idx % len(external_markers)]

                        ax.plot(
                            phi_ds,
                            eps_v_ds,
                            label=label,
                            color=color,
                            linestyle=":",
                            linewidth=2,
                            marker=marker,
                            markersize=6,
                            markerfacecolor="none",
                        )

            ax.set_xlim([-5, 365])
            ax.set_xticks(np.arange(0, 370, 45))
            ax.grid(True, which="both", linestyle="dashed", linewidth=0.5, alpha=0.6)
            ax.set_xlabel(r"$\theta\;[^\circ]$")
            ax.set_ylabel(r"$\varepsilon_{\mathrm{vol}}\; [\%]$")
            ax.set_ylim(ylim_range or [-0.10, 0.00])

        Path(self.output_dir or ".").mkdir(parents=True, exist_ok=True)

        if layout == "single":
            handles, labels = ax.get_legend_handles_labels()
            fig.subplots_adjust(right=0.75)
            fig.legend(
                handles,
                labels,
                loc="center left",
                bbox_to_anchor=(1.01, 0.5),
                frameon=False,
                fontsize="small",
            )
            plt.savefig(
                Path(self.output_dir, f"volume_strain_plot_{model_type}_all.png"),
                dpi=400,
                bbox_inches="tight",
            )
            plt.show()
            plt.close()

        elif layout == "subplots":
            fig.subplots_adjust(right=0.75, hspace=0.4)

            for load_label in load_labels:
                ax = axes[load_label]
                handles, labels = ax.get_legend_handles_labels()
                ax.legend(
                    handles,
                    labels,
                    loc="center left",
                    bbox_to_anchor=(1.01, 0.5),
                    fontsize="small",
                    frameon=False,
                )

            plot_path = Path(
                self.output_dir, f"volume_strain_plot_{model_type}_subplots.png"
            )
            plt.savefig(plot_path, dpi=400, bbox_inches="tight")
            plt.show()
            plt.close()
