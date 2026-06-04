from __future__ import annotations

from pathlib import Path

import numpy as np
import ogstools as ot
import pyvista as pv
from roughfracture_runner import jrc_dir_name, sigma_tag

CMAP_VMAG = "jet"
CMAP_P = "turbo"
CMAP_AP = "viridis"
CMAP_CT = "gray"

n_stream_seeds = 50
boundary_strip_pct = 0.03
stream_clim_pct = (2, 98)
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
            return None, None

    v_all = np.asarray(surf.point_data["v"])
    vm_all = _vmag(v_all)
    surf.point_data["v_mag"] = vm_all

    bounds = surf.bounds
    int_clim = _interior_pct_clim(vm_all, np.asarray(surf.points), bounds)
    diag = float(surf.length)

    dx = bounds[1] - bounds[0]
    x_in = bounds[0] + dx * 0.01
    y_seed = np.linspace(bounds[2], bounds[3], n_stream_seeds + 2)[1:-1]
    source = pv.PolyData(
        np.column_stack([np.full(len(y_seed), x_in), y_seed, np.zeros(len(y_seed))])
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
        return None, int_clim

    if sl is None or sl.n_points == 0:
        return None, int_clim

    sl_v = np.asarray(sl.point_data.get("v", np.zeros((sl.n_points, 3))))
    sl_vm = np.linalg.norm(sl_v, axis=1)
    sl_vm[sl_vm <= 0] = min_mag
    sl.point_data["v_mag"] = sl_vm

    return sl, int_clim


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
                    open_vals = vals[vals > 1e-7]
                    if open_vals.size > 10:
                        vmins.append(0.0)
                        vmaxs.append(float(np.percentile(open_vals, 98)))
                    elif vals.size:
                        vmins.append(0.0)
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


def render_scalar_grid(
    jrc_list: list[int],
    sigmas_mpa: list[float],
    scalar_name: str,
    cmap_name: str,
    cbar_label: str,
    base_run_dir: Path,
    zoom_factor: float = 1.4,
    cell_size: tuple = (680, 520),
) -> None:
    clim = _compute_global_clim(jrc_list, sigmas_mpa, scalar_name, base_run_dir)
    print(f"  global clim for {scalar_name}: {clim}")

    nrows, ncols = len(jrc_list), len(sigmas_mpa)
    pl = pv.Plotter(
        shape=(nrows, ncols),
        border=False,
        window_size=(ncols * cell_size[0], nrows * cell_size[1]),
    )

    for i, jrc in enumerate(jrc_list):
        for j, sigma in enumerate(sigmas_mpa):
            pl.subplot(i, j)
            pl.set_background("white")
            pl.add_text(
                f"JRC={int(jrc)},  σₙ={float(sigma):g} MPa",
                font_size=9,
                color="black",
            )

            try:
                jrc_, sigma_ = int(jrc), float(sigma)
                pvd = _find_pvd(base_run_dir, jrc_, sigma_)
                out_mesh = _load_last_output(pvd)
                in_mesh = ot.mesh.read(str(_find_input_vtu(base_run_dir, jrc_, sigma_)))
                surf = in_mesh.extract_surface(algorithm="dataset_surface")
                bg = surf.copy(deep=True)

                actual_clim = clim
                scalar_preference = "point"
                bg_opacity = 0.9

                if scalar_name == "aperture_closed":
                    scalar_preference = "cell"
                    name = _find_first(APERTURE_CANDIDATES, surf.cell_data.keys())
                    if name is None:
                        msg = "No aperture field"
                        raise KeyError(msg)
                    bg.cell_data["aperture_closed"] = np.asarray(
                        surf.cell_data[name], float
                    )

                elif scalar_name == "contact_closed":
                    scalar_preference = "cell"
                    name = _find_first(CONTACT_CANDIDATES, surf.cell_data.keys())
                    if name is None:
                        msg = "No contact field"
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
                        msg = "No pressure field"
                        raise KeyError(msg)
                    bg["pressure"] = _clean_pressure(bg.point_data[pname])
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
                        actual_clim = (
                            float(np.nanpercentile(p_int, 2)),
                            float(np.nanpercentile(p_int, 98)),
                        )

                show_bar = i == nrows - 1 and j == ncols // 2
                pl.add_mesh(
                    bg,
                    scalars=scalar_name,
                    cmap=cmap_name,
                    clim=actual_clim,
                    opacity=bg_opacity,
                    preference=scalar_preference,
                    show_edges=False,
                    show_scalar_bar=show_bar,
                    scalar_bar_args={
                        "title": cbar_label,
                        "vertical": False,
                        "position_x": 0.1,
                        "position_y": 0.03,
                        "width": 0.8,
                        "height": 0.16,
                        "title_font_size": 20,
                        "label_font_size": 16,
                        "color": "black",
                    },
                )

                sl_lines, sl_clim = _build_streamlines(out_mesh)
                if sl_lines is not None:
                    pl.add_mesh(
                        sl_lines,
                        scalars="v_mag",
                        cmap=CMAP_VMAG,
                        clim=sl_clim,
                        render_lines_as_tubes=True,
                        line_width=3.0,
                        show_scalar_bar=False,
                        lighting=False,
                    )

                pl.enable_parallel_projection()
                b = bg.bounds
                cx_ = (b[0] + b[1]) / 2
                cy_ = (b[2] + b[3]) / 2
                cz_ = (b[4] + b[5]) / 2
                diag = float(bg.length)
                pl.camera.position = (
                    cx_ + diag * 0.6,
                    cy_ + diag * 0.8,
                    cz_ + diag * 0.5,
                )
                pl.camera.focal_point = (cx_, cy_, cz_)
                pl.camera.up = (0, 0, 1)
                pl.camera.zoom(float(zoom_factor))

            except (FileNotFoundError, KeyError) as exc:
                pl.add_text(f"missing\n{exc}", font_size=8, color="red")

    pl.show()


def run_all_plots(
    *,
    out_dir: Path,
    jrc_list: list[int],
    sigmas_mpa: list[float],
) -> None:
    base_run_dir = Path(out_dir) / "runs"

    configs = [
        ("pressure", CMAP_P, r"p / Pa"),
        ("aperture_closed", CMAP_AP, r"w / m"),
        ("contact_closed", CMAP_CT, r"Ic / -"),
    ]

    for scalar_name, cmap, cbar_label in configs:
        print(f"\nPlotting {scalar_name} …")
        try:
            render_scalar_grid(
                jrc_list=jrc_list,
                sigmas_mpa=sigmas_mpa,
                scalar_name=scalar_name,
                cmap_name=cmap,
                cbar_label=cbar_label,
                base_run_dir=base_run_dir,
            )
        except Exception as exc:
            print(f"  FAILED {scalar_name}: {exc}")


def run_roughfracture_postprocessing(
    *,
    out_dir: Path,
    jrc_list: list[int],
    sigmas_mpa: list[float],
) -> None:
    run_all_plots(
        out_dir=out_dir,
        jrc_list=jrc_list,
        sigmas_mpa=sigmas_mpa,
    )
