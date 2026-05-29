from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
import pyvista as pv
from IPython.display import Image as _Img
from IPython.display import display as _display
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from roughfracture_runner import sigma_tag
from vtkmodules.vtkRenderingCore import vtkCoordinate

APERTURE_CANDIDATES = ("aperture_closed", "aperture", "w")
CONTACT_CANDIDATES = ("contact_closed", "is_contact", "contact_flag")


def _find_first(candidates, names):
    for n in candidates:
        if n in names:
            return n
    return None


def _show(fig: plt.Figure, save_to: Path | None) -> None:
    if save_to is not None:
        fig.savefig(save_to, dpi=200, bbox_inches="tight")
        print(f"saved: {save_to}")
    if "ipykernel" in sys.modules:
        _display(fig)
        plt.close(fig)
    else:
        plt.show()


def _add_boundary_overlays(pl: pv.Plotter, mesh: pv.DataSet) -> list:
    b = mesh.bounds
    tol_x = (b[1] - b[0]) * 0.015
    tol_y = (b[3] - b[2]) * 0.015
    cx = (b[0] + b[1]) / 2
    cy = (b[2] + b[3]) / 2
    z_label = b[5]

    edges = mesh.extract_feature_edges(
        boundary_edges=True,
        feature_edges=False,
        manifold_edges=False,
        non_manifold_edges=False,
    )
    pts = np.asarray(edges.points)

    sides = [
        (
            np.abs(pts[:, 0] - b[0]) < tol_x,
            "#d62728",
            "Inlet ($p_\\mathrm{in}$)",
            [b[0], cy, z_label],
        ),
        (
            np.abs(pts[:, 0] - b[1]) < tol_x,
            "#1f77b4",
            "Outlet ($p_\\mathrm{out}$)",
            [b[1], cy, z_label],
        ),
        (np.abs(pts[:, 1] - b[2]) < tol_y, "#7f7f7f", "No-flow", [cx, b[2], z_label]),
        (np.abs(pts[:, 1] - b[3]) < tol_y, "#7f7f7f", "No-flow", [cx, b[3], z_label]),
    ]

    label_specs = []
    for mask, color, label, lpt in sides:
        idx = np.where(mask)[0]
        if idx.size < 2:
            continue
        sub = pts[idx]
        axis = 0 if np.ptp(sub[:, 0]) > np.ptp(sub[:, 1]) else 1
        sub = sub[np.argsort(sub[:, axis])]
        pl.add_mesh(
            pv.Spline(sub, n_points=max(60, len(sub))),
            color=color,
            line_width=6,
            render_lines_as_tubes=True,
        )
        label_specs.append((lpt, label, color))
    return label_specs


def _project_labels(
    label_specs: list, renderer, window_size: tuple, scale: int
) -> list:
    coord = vtkCoordinate()
    coord.SetCoordinateSystemToWorld()
    win_h = window_size[1]
    projected = []
    for world_pt, text, color in label_specs:
        coord.SetValue(*world_pt)
        px, py = coord.GetComputedDoubleDisplayValue(renderer)
        projected.append((px * scale, (win_h - py) * scale, text, color))
    return projected


def _stamp_boundary_labels(png_path: Path, projected: list) -> None:

    img = plt.imread(str(png_path))
    h, w = img.shape[:2]

    fig = Figure(figsize=(w / 150, h / 150), dpi=150)
    FigureCanvasAgg(fig)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(img)
    ax.axis("off")

    fontsize = max(10, w // 110)
    for ix, iy, text, color in projected:
        ax.text(
            ix,
            iy,
            text,
            color=color,
            fontsize=fontsize,
            fontweight="bold",
            ha="center",
            va="center",
            bbox={
                "boxstyle": "round,pad=0.35",
                "facecolor": "white",
                "alpha": 0.85,
                "edgecolor": color,
                "linewidth": 1.5,
            },
        )

    fig.savefig(str(png_path), dpi=150, bbox_inches=None, pad_inches=0)


def _render_surface_png(
    vtu_path: Path,
    png_path: Path,
    *,
    z_scale: float,
    scalar: str,
    clim: tuple[float, float],
    cmap: str,
    window_size: tuple,
    show_boundaries: bool = False,
) -> None:
    mesh = pv.read(str(vtu_path))

    pts = mesh.points.copy()
    pts[:, 2] *= z_scale
    mesh.points = pts

    if scalar in mesh.cell_data:
        mesh = mesh.cell_data_to_point_data()

    pl = pv.Plotter(off_screen=True, window_size=window_size)
    pl.set_background("white")
    pl.add_mesh(
        mesh,
        scalars=scalar,
        cmap=cmap,
        clim=clim,
        show_edges=False,
        show_scalar_bar=False,
    )
    label_specs = _add_boundary_overlays(pl, mesh) if show_boundaries else []
    pl.add_text(
        f"z x{int(z_scale)}", position="lower_right", font_size=14, color="dimgray"
    )
    pl.add_axes(
        line_width=3,
        xlabel="x",
        ylabel="y",
        zlabel="z",
        x_color="#d62728",
        y_color="#2ca02c",
        z_color="#1f77b4",
    )
    pl.enable_parallel_projection()
    bounds = mesh.bounds
    cx = (bounds[0] + bounds[1]) / 2
    cy = (bounds[2] + bounds[3]) / 2
    cz = (bounds[4] + bounds[5]) / 2
    diag = float(mesh.length)
    pl.camera.position = (cx + diag * 0.6, cy + diag * 0.8, cz + diag * 0.5)
    pl.camera.focal_point = (cx, cy, cz)
    pl.camera.up = (0, 0, 1)
    pl.camera.zoom(1.2)
    pl.screenshot(str(png_path), scale=2)
    projected = _project_labels(label_specs, pl.renderer, window_size, scale=2)
    pl.close()
    if projected:
        _stamp_boundary_labels(png_path, projected)


def plot_fracture_surface_3d(
    vtu_path: Path,
    *,
    z_scale: float = 2.0,
    scalar: str = "aperture_closed",
    cmap: str = "viridis",
    window_size: tuple = (1600, 900),
    show_boundaries: bool = True,
    save_to: Path | None = None,
) -> None:
    mesh = pv.read(str(Path(vtu_path)))
    vals = (
        np.asarray(mesh.cell_data[scalar], float)
        if scalar in mesh.cell_data
        else np.zeros(mesh.n_cells)
    )
    clim = (float(vals.min()), float(vals.max()))

    if save_to is None:
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as _f:
            tmp = Path(_f.name)
    else:
        tmp = save_to
    _render_surface_png(
        Path(vtu_path),
        tmp,
        z_scale=z_scale,
        scalar=scalar,
        clim=clim,
        cmap=cmap,
        window_size=window_size,
        show_boundaries=show_boundaries,
    )
    if "ipykernel" in sys.modules:
        _display(_Img(str(tmp)))


def plot_aperture_grid(
    *,
    jrc_list: list[int],
    sigmas_mpa: list[float],
    mesh_dir: Path,
    fontsize: int = 13,
    save_to: Path | None = None,
) -> None:
    mesh_dir = Path(mesh_dir)
    nrows, ncols = len(jrc_list), len(sigmas_mpa)

    meshes: dict[tuple, pv.DataSet | None] = {}
    all_vals: list[float] = []

    for jrc in jrc_list:
        for sigma in sigmas_mpa:
            s = sigma_tag(float(sigma))
            vtu = mesh_dir / f"joint_JRC{int(jrc)}_sigma_{s}MPa.vtu"
            if not vtu.exists():
                meshes[(int(jrc), float(sigma))] = None
                continue
            m = pv.read(str(vtu))
            meshes[(int(jrc), float(sigma))] = m
            ap = _find_first(APERTURE_CANDIDATES, m.cell_data.keys())
            if ap:
                v = np.asarray(m.cell_data[ap], float)
                all_vals.extend(v[np.isfinite(v)].tolist())

    vmin = float(np.percentile(all_vals, 2))
    vmax = float(np.percentile(all_vals, 98))

    fig = plt.figure(figsize=(4.2 * ncols, 3.6 * nrows + 0.6))
    gs = fig.add_gridspec(
        nrows + 1, ncols, height_ratios=[1.0] * nrows + [0.04], hspace=0.12, wspace=0.06
    )

    for i, jrc in enumerate(jrc_list):
        for j, sigma in enumerate(sigmas_mpa):
            ax = fig.add_subplot(gs[i, j])
            key = (int(jrc), float(sigma))
            m = meshes.get(key)

            ax.set_title(
                rf"JRC = {int(jrc)},  $\sigma_n = {float(sigma):g}$ MPa",
                fontsize=fontsize,
            )
            ax.set_aspect("equal")
            ax.axis("off")

            if m is None:
                ax.text(
                    0.5,
                    0.5,
                    "missing",
                    ha="center",
                    va="center",
                    transform=ax.transAxes,
                )
                continue

            # UnstructuredGrid cells: [3, v0, v1, v2, 3, v0, v1, v2, ...]
            faces = m.cells.reshape(-1, 4)[:, 1:]
            pts = m.points
            triang = mtri.Triangulation(pts[:, 0], pts[:, 1], faces)

            mp = m.cell_data_to_point_data()
            ap_name = _find_first(APERTURE_CANDIDATES, mp.point_data.keys())
            apt = np.asarray(mp.point_data[ap_name], float)

            ax.tripcolor(
                triang, apt, cmap="viridis", vmin=vmin, vmax=vmax, shading="gouraud"
            )

            if i == nrows - 1 and j == 0:
                kw = {
                    "xycoords": "axes fraction",
                    "textcoords": "axes fraction",
                    "arrowprops": {
                        "arrowstyle": "-|>",
                        "color": "k",
                        "lw": 2.0,
                        "mutation_scale": 14,
                    },
                }
                ax.annotate("", xy=(0.22, 0.07), xytext=(0.05, 0.07), **kw)
                ax.annotate("", xy=(0.05, 0.26), xytext=(0.05, 0.07), **kw)
                ax.text(
                    0.24,
                    0.07,
                    "x",
                    transform=ax.transAxes,
                    fontsize=fontsize + 1,
                    fontweight="bold",
                    ha="left",
                    va="center",
                )
                ax.text(
                    0.05,
                    0.28,
                    "y",
                    transform=ax.transAxes,
                    fontsize=fontsize + 1,
                    fontweight="bold",
                    ha="center",
                    va="bottom",
                )

    cax = fig.add_subplot(gs[-1, :])
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    sm = mpl.cm.ScalarMappable(norm=norm, cmap="viridis")
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cbar.set_label(r"$w$ / m", fontsize=fontsize)
    cbar.ax.tick_params(labelsize=fontsize - 2)

    _show(fig, save_to)


def plot_fracture_geometry_summary(
    *,
    jrc_list: list[int],
    sigmas_mpa: list[float],
    mesh_dir: Path,
    fontsize: int = 13,
    figsize: tuple = (10, 4.0),
    save_to: Path | None = None,
) -> None:
    mesh_dir = Path(mesh_dir)
    colors = plt.cm.tab10.colors

    data: dict[int, dict] = {}
    for jrc in jrc_list:
        data[int(jrc)] = {"sigma": [], "mean_aperture": [], "contact_pct": []}
        for sigma in sigmas_mpa:
            s = sigma_tag(float(sigma))
            vtu = mesh_dir / f"joint_JRC{int(jrc)}_sigma_{s}MPa.vtu"
            if not vtu.exists():
                print(f"  skip: {vtu.name}")
                continue
            mesh = pv.read(str(vtu))

            ap_name = _find_first(APERTURE_CANDIDATES, mesh.cell_data.keys())
            if ap_name is not None:
                w = np.asarray(mesh.cell_data[ap_name], float)
                w = w[np.isfinite(w)]
                mean_w = float(np.mean(w)) if w.size else np.nan
            else:
                mean_w = np.nan

            ct_name = _find_first(CONTACT_CANDIDATES, mesh.cell_data.keys())
            contact_pct = (
                float(np.mean(np.asarray(mesh.cell_data[ct_name], float) > 0.5) * 100.0)
                if ct_name is not None
                else np.nan
            )

            data[int(jrc)]["sigma"].append(float(sigma))
            data[int(jrc)]["mean_aperture"].append(mean_w)
            data[int(jrc)]["contact_pct"].append(contact_pct)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

    for idx, jrc in enumerate(jrc_list):
        d = data[int(jrc)]
        col = colors[idx % len(colors)]
        ax1.plot(
            d["sigma"],
            d["mean_aperture"],
            "o-",
            color=col,
            label=f"JRC {int(jrc)}",
            linewidth=1.8,
        )
        ax2.plot(
            d["sigma"],
            d["contact_pct"],
            "o-",
            color=col,
            label=f"JRC {int(jrc)}",
            linewidth=1.8,
        )

    for ax, ylabel in zip(
        (ax1, ax2), (r"$\bar{w}$ / m", "contact area / %"), strict=False
    ):
        ax.set_xlabel(r"$\sigma_n$ / MPa", fontsize=fontsize)
        ax.set_ylabel(ylabel, fontsize=fontsize)
        ax.legend(fontsize=fontsize - 1, framealpha=0.4)
        ax.tick_params(labelsize=fontsize - 2)
        ax.grid(True, linestyle="--", linewidth=0.6, alpha=0.5)
        ax.spines[["top", "right"]].set_visible(False)

    fig.tight_layout()
    _show(fig, save_to)
