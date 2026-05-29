from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import ogstools as ot
import pyvista as pv
from IPython.display import display as _display
from roughfracture_runner import jrc_dir_name, sigma_tag
from scipy.spatial import cKDTree

CMAP_VMAG = "jet"
CMAP_P = "turbo"
CMAP_AP = "viridis"
CMAP_CT = "gray"

n_stream_seeds = 50
stream_tube_r_frac = 0.001
stream_tube_min_radius = 2.0e-5  # [m]
stream_tube_n_sides = 8
stream_line_width = 1.5
boundary_strip_pct = 0.03
stream_clim_pct = (2, 98)
streamline_tube_opacity = 1.0

opacity = 0.5
min_mag = 1e-14

PRESSURE_CANDIDATES = ("pressure", "LiquidFlow_pressure", "p")
APERTURE_CANDIDATES = ("aperture_closed", "aperture", "w")
CONTACT_CANDIDATES = ("contact_closed", "is_contact", "contact_flag")


def _vmag(v: np.ndarray) -> np.ndarray:
    v = np.asarray(v)
    m = np.linalg.norm(v, axis=1)
    m[m <= 0] = min_mag
    return m


def _find_first(candidates, names):
    for n in candidates:
        if n in names:
            return n
    return None


def _finite_vals(a: np.ndarray) -> np.ndarray:
    a = np.asarray(a, float).ravel()
    return a[np.isfinite(a)]


def _clean_pressure(pvals) -> np.ndarray:
    p = np.asarray(pvals, float)
    p[p <= 0.0] = np.nan
    return p


def _find_pvd(base_run_dir: Path, jrc: int, sigma_mpa: float) -> Path:
    out_dir = base_run_dir / jrc_dir_name(jrc) / f"sigma_{sigma_tag(sigma_mpa)}" / "out"
    pvds = sorted(out_dir.glob("*.pvd"))
    if not pvds:
        msg = f"No .pvd found in {out_dir}"
        raise FileNotFoundError(msg)
    return pvds[0]


def _find_input_vtu(base_run_dir: Path, jrc: int, sigma_mpa: float) -> Path:
    mesh_dir = (
        base_run_dir / jrc_dir_name(jrc) / f"sigma_{sigma_tag(sigma_mpa)}" / "mesh"
    )
    skip = {"left.vtu", "right.vtu", "top.vtu", "bottom.vtu", "boundaries.vtu"}
    for vtu in sorted(mesh_dir.glob("*.vtu")):
        if vtu.name not in skip:
            return vtu
    msg = f"No domain mesh in {mesh_dir}"
    raise FileNotFoundError(msg)


def _load_last_output(pvd: Path) -> pv.DataSet:
    series = ot.MeshSeries(str(pvd))
    mesh = series[len(series) - 1]
    if len(mesh.cell_data.keys()) > 0:
        mesh = mesh.cell_data_to_point_data(pass_cell_data=False)
    return mesh


def _interior_pct_clim(
    v_mag: np.ndarray, points: np.ndarray, bounds
) -> tuple[float, float]:
    dx = bounds[1] - bounds[0]
    dy = bounds[3] - bounds[2]
    sx, sy = dx * boundary_strip_pct, dy * boundary_strip_pct
    mask = (
        (points[:, 0] > bounds[0] + sx)
        & (points[:, 0] < bounds[1] - sx)
        & (points[:, 1] > bounds[2] + sy)
        & (points[:, 1] < bounds[3] - sy)
    )
    interior = v_mag[mask]
    base = interior if interior.size > 20 else v_mag
    return (
        float(np.percentile(base, stream_clim_pct[0])),
        float(np.percentile(base, stream_clim_pct[1])),
    )


def _build_streamlines(output_mesh: pv.DataSet):
    """Return (tubes, interior_clim, used_tubes) from output velocity field."""
    surf = (
        output_mesh.copy(deep=False)
        if isinstance(output_mesh, pv.PolyData)
        else output_mesh.extract_surface(algorithm="dataset_surface")
    )

    if "v" not in surf.point_data:
        if "v" in surf.cell_data:
            surf = surf.cell_data_to_point_data()
        elif "v" in surf.array_names:
            surf.point_data["v"] = surf["v"]
        else:
            return None, None, False

    v_all = np.asarray(surf.point_data["v"])
    vm_all = _vmag(v_all)
    surf.point_data["v_mag"] = vm_all

    bounds = surf.bounds
    int_clim = _interior_pct_clim(vm_all, np.asarray(surf.points), bounds)
    diag = float(surf.length)
    tube_r = max(diag * stream_tube_r_frac, stream_tube_min_radius)

    dx = bounds[1] - bounds[0]
    x_in = bounds[0] + dx * 0.01
    y_seed = np.linspace(bounds[2], bounds[3], n_stream_seeds + 2)[1:-1]
    source = pv.PolyData(
        np.column_stack(
            [
                np.full(len(y_seed), x_in),
                y_seed,
                np.zeros(len(y_seed)),
            ]
        )
    )

    try:
        sl = surf.streamlines_from_source(
            source,
            vectors="v",
            surface_streamlines=True,
            max_steps=5000,
            initial_step_length=diag * 0.004,
            terminal_speed=min_mag * 10,
            integration_direction="forward",
        )
    except Exception as exc:
        print(f"  [streamlines] failed: {exc}")
        return None, int_clim, False

    if sl is None or sl.n_points == 0:
        return None, int_clim, False

    sl_v = np.asarray(sl.point_data.get("v", np.zeros((sl.n_points, 3))))
    sl_vm = np.linalg.norm(sl_v, axis=1)
    sl_vm[sl_vm <= 0] = min_mag
    sl.point_data["v_mag"] = sl_vm

    try:
        tubes = sl.tube(radius=tube_r, n_sides=stream_tube_n_sides, capping=False)
        if tubes is None or tubes.n_points == 0:
            msg = "tube() returned empty geometry"
            raise RuntimeError(msg)
        _, idx = cKDTree(sl.points).query(tubes.points)
        tubes.point_data["v_mag"] = sl_vm[idx]
        used_tubes = True
    except Exception as exc:
        print(f"  [streamlines] tube failed ({exc}); using polylines")
        tubes = sl
        used_tubes = False

    return tubes, int_clim, used_tubes


def render_case_to_png(
    jrc: int,
    sigma: float,
    png_file: Path,
    scalar_name: str,
    clim: tuple,
    cmap_name: str,
    base_run_dir: Path,
    window_size: tuple = (2600, 1800),
    scale: int = 2,
    zoom_factor: float = 1.0,
) -> None:
    jrc, sigma = int(jrc), float(sigma)
    pvd = _find_pvd(base_run_dir, jrc, sigma)
    out_mesh = _load_last_output(pvd)
    in_mesh = ot.mesh.read(str(_find_input_vtu(base_run_dir, jrc, sigma)))
    surf = in_mesh.extract_surface(algorithm="dataset_surface")
    bg = surf.copy(deep=True)

    scalar_preference = "point"
    bg_opacity = float(opacity)

    if scalar_name in ("aperture_closed", "contact_closed"):
        scalar_preference = "cell"
        if scalar_name == "aperture_closed":
            name = _find_first(APERTURE_CANDIDATES, surf.cell_data.keys())
            if name is None:
                msg = (
                    f"No aperture field in {_find_input_vtu(base_run_dir, jrc, sigma)}"
                )
                raise KeyError(msg)
            bg.cell_data["aperture_closed"] = np.asarray(surf.cell_data[name], float)
        else:
            name = _find_first(CONTACT_CANDIDATES, surf.cell_data.keys())
            if name is None:
                msg = f"No contact field in {_find_input_vtu(base_run_dir, jrc, sigma)}"
                raise KeyError(msg)
            bg.cell_data["contact_closed"] = (
                np.asarray(surf.cell_data[name], float) > 0.5
            ).astype(float)
            bg_opacity = 1.0

    elif scalar_name == "pressure":
        sampled = out_mesh.sample(surf)
        for k in list(sampled.point_data.keys()):
            arr = np.asarray(sampled.point_data[k])
            if arr.shape[0] == bg.n_points:
                bg.point_data[k] = arr

        pname = _find_first(PRESSURE_CANDIDATES, bg.point_data.keys())
        if pname is None:
            msg = "No pressure field in output"
            raise KeyError(msg)
        bg["pressure"] = _clean_pressure(bg.point_data[pname])
        # use per-case interior percentile clim (pressure range differs across cases)
        p_vals = np.asarray(bg["pressure"], float)
        b_ = bg.bounds
        pts_ = np.asarray(bg.points)
        sx_ = (b_[1] - b_[0]) * boundary_strip_pct
        sy_ = (b_[3] - b_[2]) * boundary_strip_pct
        mask_ = (
            (pts_[:, 0] > b_[0] + sx_)
            & (pts_[:, 0] < b_[1] - sx_)
            & (pts_[:, 1] > b_[2] + sy_)
            & (pts_[:, 1] < b_[3] - sy_)
        )
        p_int = p_vals[mask_]
        p_int = p_int[np.isfinite(p_int)]
        if p_int.size > 10:
            clim = (
                float(np.nanpercentile(p_int, 2)),
                float(np.nanpercentile(p_int, 98)),
            )
    else:
        msg = f"Unsupported scalar_name='{scalar_name}'"
        raise ValueError(msg)

    sl_tubes, sl_clim, sl_used_tubes = _build_streamlines(out_mesh)

    p = pv.Plotter(off_screen=True, window_size=window_size)
    p.set_background("white")
    p.add_mesh(
        bg,
        scalars=scalar_name,
        cmap=cmap_name,
        clim=clim,
        opacity=bg_opacity,
        preference=scalar_preference,
        show_edges=False,
        show_scalar_bar=False,
    )

    if sl_tubes is not None:
        sl_kw = {
            "scalars": "v_mag",
            "cmap": CMAP_VMAG,
            "clim": sl_clim,
            "opacity": streamline_tube_opacity,
            "show_scalar_bar": False,
            "lighting": False,
        }
        if not sl_used_tubes:
            sl_kw["render_lines_as_tubes"] = False
            sl_kw["line_width"] = float(stream_line_width)
        p.add_mesh(sl_tubes, **sl_kw)
    p.add_axes(
        line_width=3,
        xlabel="x",
        ylabel="y",
        zlabel="z",
        x_color="#d62728",
        y_color="#2ca02c",
        z_color="#1f77b4",
    )
    p.enable_parallel_projection()
    bounds = bg.bounds
    cx = (bounds[0] + bounds[1]) / 2
    cy = (bounds[2] + bounds[3]) / 2
    cz = (bounds[4] + bounds[5]) / 2
    diag = float(bg.length)
    p.camera.position = (cx + diag * 0.6, cy + diag * 0.8, cz + diag * 0.5)
    p.camera.focal_point = (cx, cy, cz)
    p.camera.up = (0, 0, 1)
    p.camera.zoom(float(zoom_factor))
    p.screenshot(str(png_file), scale=scale)
    p.close()


def _compute_global_clim(
    jrc_list: list[int],
    sigmas_mpa: list[float],
    scalar_name: str,
    base_run_dir: Path,
) -> tuple[float, float]:
    if scalar_name == "contact_closed":
        return (0.0, 1.0)

    vmins: list[float] = []
    vmaxs: list[float] = []

    for jrc in jrc_list:
        for sigma in sigmas_mpa:
            try:
                surf = ot.mesh.read(
                    str(_find_input_vtu(base_run_dir, int(jrc), float(sigma)))
                ).extract_surface(algorithm="dataset_surface")

                if scalar_name == "aperture_closed":
                    name = _find_first(APERTURE_CANDIDATES, surf.cell_data.keys())
                    if name is None:
                        continue
                    vals = _finite_vals(np.asarray(surf.cell_data[name], float))
                    if vals.size:
                        vmins.append(float(np.min(vals)))
                        vmaxs.append(float(np.max(vals)))

                elif scalar_name == "pressure":
                    pvd = _find_pvd(base_run_dir, int(jrc), float(sigma))
                    out_mesh = _load_last_output(pvd)
                    sampled = out_mesh.sample(surf)
                    pname = _find_first(PRESSURE_CANDIDATES, sampled.array_names)
                    if pname is None:
                        continue
                    p_all = _clean_pressure(np.asarray(sampled[pname], float))
                    b_ = surf.bounds
                    pts_ = np.asarray(surf.points)
                    sx_ = (b_[1] - b_[0]) * boundary_strip_pct
                    sy_ = (b_[3] - b_[2]) * boundary_strip_pct
                    mask_ = (
                        (pts_[:, 0] > b_[0] + sx_)
                        & (pts_[:, 0] < b_[1] - sx_)
                        & (pts_[:, 1] > b_[2] + sy_)
                        & (pts_[:, 1] < b_[3] - sy_)
                    )
                    p_int = p_all[mask_]
                    p_int = p_int[np.isfinite(p_int)]
                    if p_int.size > 10:
                        vmins.append(float(np.nanpercentile(p_int, stream_clim_pct[0])))
                        vmaxs.append(float(np.nanpercentile(p_int, stream_clim_pct[1])))

            except FileNotFoundError:
                continue

    if not vmins:
        msg = f"No valid values for '{scalar_name}' under {base_run_dir}"
        raise RuntimeError(msg)
    return (float(np.min(vmins)), float(np.max(vmaxs)))


def _render_scalar_grid(
    jrc_list: list[int],
    sigmas_mpa: list[float],
    scalar_name: str,
    cmap_name: str,
    base_run_dir: Path,
    out_img_dir: Path,
    zoom_factor: float = 1.4,
) -> tuple[dict, tuple]:
    out_dir = out_img_dir / scalar_name
    out_dir.mkdir(parents=True, exist_ok=True)

    clim = _compute_global_clim(jrc_list, sigmas_mpa, scalar_name, base_run_dir)
    print(f"  global clim for {scalar_name}: {clim}")

    png_map: dict = {}
    for jrc in jrc_list:
        for sigma in sigmas_mpa:
            s = sigma_tag(float(sigma))
            png_file = out_dir / f"JRC{int(jrc)}_sigma_{s}__{scalar_name}.png"
            try:
                render_case_to_png(
                    jrc=int(jrc),
                    sigma=float(sigma),
                    png_file=png_file,
                    scalar_name=scalar_name,
                    clim=clim,
                    cmap_name=cmap_name,
                    base_run_dir=base_run_dir,
                    window_size=(2800, 2000),
                    scale=2,
                    zoom_factor=zoom_factor,
                )
                png_map[(int(jrc), float(sigma))] = png_file
                print(f"    saved: {png_file.name}")
            except FileNotFoundError as exc:
                print(f"    skip  JRC={jrc}, s={sigma}: {exc}")
                png_map[(int(jrc), float(sigma))] = None
            except Exception as exc:
                print(f"    FAILED JRC={jrc}, s={sigma}: {exc}")
                png_map[(int(jrc), float(sigma))] = None
    return png_map, clim


def _plot_grid(
    png_map: dict,
    jrc_list: list[int],
    sigmas_mpa: list[float],
    clim: tuple,
    cmap_name: str,
    cbar_label: str,
    title_fontsize: int = 16,
    cbar_fontsize: int = 16,
    cbar_ticksize: int = 12,
    figsize: tuple | None = None,
) -> plt.Figure:
    nrows, ncols = len(jrc_list), len(sigmas_mpa)
    if figsize is None:
        figsize = (6 * ncols, 4.8 * nrows + 1.0)

    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(
        nrows + 1,
        ncols,
        height_ratios=[1.0] * nrows + [0.04],
        hspace=0.02,
        wspace=0.02,
    )

    for i, jrc in enumerate(jrc_list):
        for j, sigma in enumerate(sigmas_mpa):
            ax = fig.add_subplot(gs[i, j])
            png = png_map.get((int(jrc), float(sigma)))
            ax.set_title(
                rf"JRC={int(jrc)} | $\sigma_n={float(sigma):g}\ \mathrm{{MPa}}$",
                fontsize=title_fontsize,
            )
            if png is None:
                ax.text(
                    0.5,
                    0.5,
                    "missing",
                    ha="center",
                    va="center",
                    transform=ax.transAxes,
                )
            else:
                ax.imshow(plt.imread(png))
            ax.axis("off")

    cax = fig.add_subplot(gs[-1, :])
    norm = mpl.colors.Normalize(vmin=float(clim[0]), vmax=float(clim[1]))
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.colormaps[cmap_name])
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cbar.set_label(cbar_label, fontsize=cbar_fontsize)
    cbar.ax.tick_params(labelsize=cbar_ticksize)
    return fig


def run_all_plots(
    *,
    out_dir: Path,
    jrc_list: list[int],
    sigmas_mpa: list[float],
    fontsize: int = 16,
) -> dict[str, Path]:
    """Render JRC x sigma field grids: pressure + streamlines, aperture, contact."""
    base_run_dir = Path(out_dir) / "runs"
    plot_dir = Path(out_dir) / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    configs = [
        # (scalar_name,      cmap,      cbar_label,         fname)
        ("pressure", CMAP_P, r"$p$ / Pa", "GRID_pressure.png"),
        ("aperture_closed", CMAP_AP, r"$w$ / m", "GRID_aperture_closed.png"),
        ("contact_closed", CMAP_CT, r"$I_c$ / -", "GRID_contact.png"),
    ]

    if "ipykernel" in sys.modules:

        def _show_fig(f):
            _display(f)
            plt.close(f)

    else:

        def _show_fig(f):
            plt.show()
            plt.close(f)

    saved: dict[str, Path] = {}
    for scalar_name, cmap, cbar_label, fname in configs:
        print(f"\nPlotting {scalar_name} …")
        try:
            png_map, clim = _render_scalar_grid(
                jrc_list=jrc_list,
                sigmas_mpa=sigmas_mpa,
                scalar_name=scalar_name,
                cmap_name=cmap,
                base_run_dir=base_run_dir,
                out_img_dir=plot_dir / "pngs",
            )
            fig = _plot_grid(
                png_map=png_map,
                jrc_list=jrc_list,
                sigmas_mpa=sigmas_mpa,
                clim=clim,
                cmap_name=cmap,
                cbar_label=cbar_label,
                title_fontsize=fontsize,
                cbar_fontsize=fontsize,
                cbar_ticksize=fontsize - 2,
            )
            out_file = plot_dir / fname
            fig.savefig(out_file, dpi=300, bbox_inches="tight")
            _show_fig(fig)
            saved[scalar_name] = out_file
            print(f"  saved: {out_file}")
        except Exception as exc:
            print(f"  FAILED {scalar_name}: {exc}")

    return saved


def run_roughfracture_postprocessing(
    *,
    out_dir: Path,
    jrc_list: list[int],
    sigmas_mpa: list[float],
    fontsize: int = 16,
) -> dict:
    """Run all post-processing plots."""
    return {
        "plots": run_all_plots(
            out_dir=out_dir,
            jrc_list=jrc_list,
            sigmas_mpa=sigmas_mpa,
            fontsize=fontsize,
        ),
    }
