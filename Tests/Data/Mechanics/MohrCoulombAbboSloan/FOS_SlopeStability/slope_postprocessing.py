from __future__ import annotations

import os
import re
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import ogstools as ot
import pyvista as pv
from slope_utils import ssr_params, ssr_schedule

DEFAULT_MPL_RC = {
    "font.family": "serif",
    "font.serif": ["DejaVu Serif"],
    "font.size": 16,
    "axes.labelsize": 16,
    "axes.titlesize": 16,
    "legend.fontsize": 13,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "lines.linewidth": 1.6,
    "axes.linewidth": 0.8,
    "grid.linewidth": 0.5,
    "savefig.dpi": 300,
    "mathtext.fontset": "stix",
}

_TS_RE = re.compile(r"_ts_(\d+)_t_([0-9.eE+\-]+)")
_T_FROM_NAME = re.compile(r"_t_([0-9.eE+\-]+)\.vtu$")
_CASE_NAME_RE = re.compile(
    r"^(?P<basename>.+?)__ord-(?P<ord>lin|quad)__h-(?P<h>[^_]+)"
    r"__dtmin-(?P<dtmin>[^_]+)__tol-(?P<tol>[^_]+)"
    r"__load-(?P<load_tag>.+?)__(?P<label>.+)$"
)

FIGSIZE_1x3 = (18.0, 5.2)
DPI = 200
WSPACE = 0.35
CB_PAD = 0.05
CB_SIZE = 0.025


def configure_plots(*, mpl_rc=None, ogs_fontsize: int = 11, scale: float = 1.0) -> None:
    rc = DEFAULT_MPL_RC.copy()
    if mpl_rc:
        rc.update(mpl_rc)
    for k in (
        "font.size",
        "axes.labelsize",
        "axes.titlesize",
        "legend.fontsize",
        "xtick.labelsize",
        "ytick.labelsize",
    ):
        rc[k] = float(rc[k]) * float(scale)

    plt.rcParams.update(rc)
    ot.plot.setup.fontsize = int(float(ogs_fontsize) * float(scale))
    ot.plot.setup.tick_pad = 3
    ot.plot.setup.tick_length = 4
    ot.plot.setup.linewidth = 0.8
    ot.plot.setup.dpi = int(plt.rcParams.get("savefig.dpi", 300))


def parse_last_ok_and_first_fail(
    log_path: Path,
) -> tuple[int, float, int | None, float | None]:
    txt = log_path.read_text(errors="ignore")

    ok_steps = [int(s) for s in re.findall(r"Time step #(\d+) took", txt)]
    if not ok_steps:
        msg = "No finished steps found in ogs.log (pattern: 'Time step #N took')."
        raise RuntimeError(msg)
    last_ok = max(ok_steps)

    failed_steps = sorted(
        {
            int(s)
            for s in re.findall(r"The nonlinear solver failed in time step #(\d+)", txt)
        }
    )
    first_fail = next((s for s in failed_steps if s > last_ok), None)

    def step_start_time(step: int) -> float:
        pat = rf"Time step #{step} started\. Time: ([0-9.eE+\-]+)"
        m = re.search(pat, txt)
        if not m:
            msg = f"Cannot find start time for step #{step} in ogs.log."
            raise RuntimeError(msg)
        return float(m.group(1).rstrip("."))

    t_ok = step_start_time(last_ok)
    t_fail = step_start_time(first_fail) if first_fail is not None else None
    return last_ok, t_ok, first_fail, t_fail


def list_raw_vtus(results_dir: Path, prefix: str) -> list[tuple[int, float, Path]]:
    vtus: list[tuple[int, float, Path]] = []
    for p in results_dir.glob(f"{prefix}_ts_*_t_*.vtu"):
        m = _TS_RE.search(p.stem)
        if m:
            vtus.append((int(m.group(1)), float(m.group(2)), p))
    vtus.sort(key=lambda x: x[0])
    return vtus


def write_pvd_from_vtus(
    vtu_triplets: list[tuple[int, float, Path]], out_pvd: Path
) -> Path:
    out_pvd = Path(out_pvd)
    out_pvd.parent.mkdir(parents=True, exist_ok=True)

    vtkfile = ET.Element(
        "VTKFile", type="Collection", version="0.1", byte_order="LittleEndian"
    )
    collection = ET.SubElement(vtkfile, "Collection")

    base_dir = out_pvd.parent
    for _, tt, p in vtu_triplets:
        rel = os.path.relpath(Path(p).resolve(), base_dir)
        ET.SubElement(
            collection, "DataSet", timestep=f"{tt:.16g}", group="", part="0", file=rel
        )

    ET.ElementTree(vtkfile).write(out_pvd, encoding="utf-8", xml_declaration=True)

    first_rel = collection.findall("DataSet")[0].attrib["file"]
    missing_file = base_dir / first_rel
    if not missing_file.exists():
        msg = f"PVD points to missing file: {missing_file}"
        raise FileNotFoundError(msg)

    return out_pvd


def _savefig(fig: plt.Figure, out: Path | None) -> None:
    if out is not None:
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, bbox_inches="tight", pad_inches=0.05)
    plt.show()
    plt.close(fig)


def _case_attr(case: Any, name: str, default=None):
    return getattr(
        case, name, case.get(name, default) if isinstance(case, dict) else default
    )


def _flatten_any(a: Any) -> np.ndarray:
    if a is None:
        return np.asarray([])
    if isinstance(a, dict):
        parts = [np.asarray(v).reshape(-1) for v in a.values()]
        return np.concatenate(parts) if parts else np.asarray([])
    return np.asarray(a).reshape(-1)


def _decode_h_token(h_token: str) -> str:
    return h_token.replace("m", "-").replace("p", ".")


def _parse_case_name(name: str) -> dict[str, str] | None:
    m = _CASE_NAME_RE.match(name)
    if not m:
        return None
    g = m.groupdict()
    g["h"] = _decode_h_token(g["h"])
    g["element"] = "quadratic" if g["ord"] == "quad" else "linear"
    return g


def _read_time_stepping_type(prj_path: Path) -> str | None:
    if not prj_path.exists():
        return None
    try:
        root = ET.parse(prj_path).getroot()
        node = root.find(".//time_loop/processes/process/time_stepping/type")
        if node is None or node.text is None:
            return None
        txt = node.text.strip()
        return txt or None
    except Exception:
        return None


def _plot_contourf_compat(
    mesh: pv.UnstructuredGrid,
    fig: plt.Figure,
    ax: plt.Axes,
    var,
    *,
    cmap: str,
    vmin: float | None,
    vmax: float | None,
    title: str,
    fontsize: float,
    cb_labelsize: float,
    cb_pad: float,
    cb_size: float,
    cb_label: str | None = None,
) -> None:
    kwargs_full = {
        "cmap": cmap,
        "vmin": vmin,
        "vmax": vmax,
        "fontsize": fontsize,
        "cb_labelsize": cb_labelsize,
        "cb_pad": cb_pad,
        "cb_size": cb_size,
    }
    if cb_label is not None:
        kwargs_full["cb_label"] = cb_label
    try:
        ot.plot.contourf(mesh, var, fig, ax, **kwargs_full)
    except TypeError:
        ot.plot.contourf(mesh, var, fig, ax, cmap=cmap, vmin=vmin, vmax=vmax)

    if cb_label is not None and fig.axes:
        cax = fig.axes[-1]
        cax.set_ylabel(cb_label, fontsize=cb_labelsize)

    ax.set_title(title, fontsize=plt.rcParams["axes.titlesize"], pad=3)
    ax.set_aspect("equal")
    ax.tick_params(
        direction="in", top=True, right=True, labelsize=plt.rcParams["xtick.labelsize"]
    )


def plot_F_vs_time(
    times: np.ndarray,
    Fvals: np.ndarray,
    *,
    out_png: Path | None,
    t_ok: float | None,
    t_fail: float | None,
) -> None:
    fig, ax = plt.subplots(figsize=(6.8, 3.0), layout="constrained")
    ax.plot(times, Fvals)
    if t_ok is not None:
        ax.axvline(t_ok, linestyle="--", linewidth=1.0)
    if t_fail is not None:
        ax.axvline(t_fail, linestyle="--", linewidth=1.0)
    ax.set_xlabel(r"Pseudo-time $t$ / 1")
    ax.set_ylabel(r"Reduction factor $F(t)$ / 1")
    ax.grid(True, which="major", linestyle=":", alpha=0.7)
    ax.tick_params(direction="in", top=True, right=True)
    _savefig(fig, out_png)


def plot_point_disp_vs_time(
    ms: ot.MeshSeries, point_xyz: np.ndarray, *, out_png: Path | None
) -> None:
    m0 = ms.mesh(0)
    pts = np.asarray(m0.points)
    p = np.asarray(point_xyz, dtype=float)

    idx = int(np.argmin(np.linalg.norm(pts[:, :2] - p[:2], axis=1)))
    snap = pts[idx]

    probe = ot.MeshSeries.probe(ms, [snap])
    times = np.asarray(probe.timevalues)

    disp = np.asarray(ot.variables.displacement.transform(probe))[:, 0, :]
    if disp.ndim == 1:
        disp = disp[:, None]
    if disp.shape[1] == 1:
        disp = np.column_stack([disp[:, 0], np.zeros(len(disp)), np.zeros(len(disp))])
    if disp.shape[1] == 2:
        disp = np.column_stack([disp[:, 0], disp[:, 1], np.zeros(len(disp))])

    Umag = np.linalg.norm(disp[:, :2], axis=1)

    fig, ax = plt.subplots(figsize=(6.8, 3.0), layout="constrained")
    ax.plot(times, disp[:, 0], label=r"$u_x$")
    ax.plot(times, disp[:, 1], label=r"$u_y$")
    ax.plot(times, Umag, label=r"$\|\mathbf{u}\|$")

    ax.set_xlabel(r"Pseudo-time $t$ / 1")
    ax.set_ylabel(r"Displacement / m")
    ax.grid(True, which="major", linestyle=":", alpha=0.7)
    ax.legend(frameon=True, framealpha=0.95, edgecolor="0.3")
    ax.tick_params(direction="in", top=True, right=True)
    _savefig(fig, out_png)


def _select_stress_mesh(
    mesh: pv.UnstructuredGrid, *, prefer_ip: bool
) -> tuple[pv.UnstructuredGrid, str]:
    pd = mesh.point_data or {}
    cd = mesh.cell_data or {}
    has_ip = ("sigma_ip" in pd) or ("sigma_ip" in cd)
    if prefer_ip and has_ip:
        try:
            return mesh.to_ip_mesh(), "sigma_ip"
        except Exception:
            return mesh, "sigma"
    return mesh, "sigma"


def find_vtu_at_time(
    vtus: list[tuple[int, float, Path]], t_target: float
) -> tuple[int, float, Path]:
    if not vtus:
        msg = "Empty VTU list."
        raise ValueError(msg)
    tt = float(t_target)
    return min(vtus, key=lambda x: abs(x[1] - tt))


def displacement_norm_relative_mm(
    mesh_end: pv.UnstructuredGrid, mesh_t0: pv.UnstructuredGrid
) -> np.ndarray:
    u = ot.variables.displacement
    u_end = np.asarray(u.transform(mesh_end))
    u_0 = np.asarray(u.transform(mesh_t0))

    u_end = np.atleast_2d(u_end)
    u_0 = np.atleast_2d(u_0)

    du = u_end - u_0
    if du.shape[1] >= 2:
        rel = np.linalg.norm(du[:, :2], axis=1)
    else:
        rel = np.abs(du[:, 0])
    return 1000.0 * rel  # mm


def _time_from_vtu_path(vtu: Path) -> float | None:
    m = _T_FROM_NAME.search(vtu.name)
    return float(m.group(1)) if m else None


def plot_baseline_stress_row_t0(
    *,
    post_dir: Path,
    vtu_t0: Path,
    prefer_ip_stress: bool = True,
) -> None:
    post_dir.mkdir(parents=True, exist_ok=True)

    t0_val = _time_from_vtu_path(vtu_t0)
    t0_txt = f"{t0_val:.5g}" if t0_val is not None else "t_0"

    mesh0 = ot.MeshSeries(str(vtu_t0))[-1]
    s_mesh, data_name = _select_stress_mesh(mesh0, prefer_ip=prefer_ip_stress)

    p_mean = ot.variables.stress.tensor_mean.replace(
        data_name=data_name, output_unit="kPa"
    )
    s = ot.variables.stress.replace(data_name=data_name, output_unit="kPa")

    def _principal_min_index(s_var, mesh: pv.UnstructuredGrid) -> int:
        l0 = np.asarray(s_var.eigenvalues[0].transform(mesh)).reshape(-1)
        l1 = np.asarray(s_var.eigenvalues[1].transform(mesh)).reshape(-1)
        l2 = np.asarray(s_var.eigenvalues[2].transform(mesh)).reshape(-1)

        lam = np.vstack([l0, l1, l2])  # (3, n)
        argmin = np.nanargmin(lam, axis=0)  # per location: 0/1/2

        counts = np.bincount(argmin, minlength=3)
        return int(np.argmax(counts))

    s = ot.variables.stress.replace(data_name=data_name, output_unit="kPa")
    idx_min = _principal_min_index(s, s_mesh)
    sig_min_var = s.eigenvalues[idx_min]

    strain_tr = ot.variables.strain.trace

    fig, axs = plt.subplots(1, 3, figsize=FIGSIZE_1x3, dpi=DPI, constrained_layout=True)
    fig.suptitle(rf"$t={t0_txt}$", y=1.02)

    _plot_contourf_compat(
        s_mesh,
        fig,
        axs[0],
        p_mean,
        cmap="coolwarm",
        vmin=None,
        vmax=None,
        title=r"$\pi$",
        fontsize=plt.rcParams["xtick.labelsize"],
        cb_labelsize=plt.rcParams["xtick.labelsize"],
        cb_pad=CB_PAD,
        cb_size=CB_SIZE,
    )

    _plot_contourf_compat(
        s_mesh,
        fig,
        axs[1],
        sig_min_var,
        cmap="coolwarm",
        vmin=None,
        vmax=None,
        title=r"$\sigma_{\min} (\lambda_0)$",
        fontsize=plt.rcParams["xtick.labelsize"],
        cb_labelsize=plt.rcParams["xtick.labelsize"],
        cb_pad=CB_PAD,
        cb_size=CB_SIZE,
    )

    _plot_contourf_compat(
        mesh0,
        fig,
        axs[2],
        strain_tr,
        cmap="RdBu_r",
        vmin=None,
        vmax=None,
        title=r"$\mathrm{tr}(\boldsymbol{\varepsilon})$",
        fontsize=plt.rcParams["xtick.labelsize"],
        cb_labelsize=plt.rcParams["xtick.labelsize"],
        cb_pad=CB_PAD,
        cb_size=CB_SIZE,
    )

    _savefig(fig, post_dir / "baseline_t0_p_sigmin_trEps.png")


def plot_end_fields_row_tEnd(
    *,
    post_dir: Path,
    vtu_t0: Path,
    vtu_end: Path,
    ranges: dict[str, tuple[float | None, float | None]] | None,
    epsp_candidates: tuple[str, ...] = (
        "EquivalentPlasticStrain",
        "eps_p",
        "epsp",
        "equivalent_plastic_strain",
    ),
) -> None:
    post_dir.mkdir(parents=True, exist_ok=True)
    ranges = ranges or {}

    tE_val = _time_from_vtu_path(vtu_end)
    tE_txt = f"{tE_val:.5g}" if tE_val is not None else r"t_{\mathrm{end}}"

    mesh0 = ot.MeshSeries(str(vtu_t0))[-1]
    meshE = ot.MeshSeries(str(vtu_end))[-1]

    u_rel = displacement_norm_relative_mm(meshE, mesh0)
    meshE.point_data["u_rel_norm_mm"] = u_rel

    strain_tr = ot.variables.strain.trace

    pd = meshE.point_data or {}
    cd = meshE.cell_data or {}
    epsp_name = next((k for k in epsp_candidates if (k in pd) or (k in cd)), None)

    def get_range(key: str):
        return ranges.get(key, (None, None))

    fig, axs = plt.subplots(1, 3, figsize=FIGSIZE_1x3, dpi=DPI, constrained_layout=True)
    fig.suptitle(rf"$t={tE_txt}$", y=1.02)

    vmin, vmax = get_range("u_rel")
    _plot_contourf_compat(
        meshE,
        fig,
        axs[0],
        "u_rel_norm_mm",
        cmap="viridis",
        vmin=vmin,
        vmax=vmax,
        title=r"$\|\mathbf{u}-\mathbf{u}_0\|$",
        fontsize=plt.rcParams["xtick.labelsize"],
        cb_labelsize=plt.rcParams["xtick.labelsize"],
        cb_pad=CB_PAD,
        cb_size=CB_SIZE,
        cb_label="Relative displacement / mm",
    )

    if epsp_name:
        vmin, vmax = get_range("epsp")
        _plot_contourf_compat(
            meshE,
            fig,
            axs[1],
            epsp_name,
            cmap="magma",
            vmin=vmin,
            vmax=vmax,
            title=r"$\bar{\varepsilon}^{\,p}$",
            fontsize=plt.rcParams["xtick.labelsize"],
            cb_labelsize=plt.rcParams["xtick.labelsize"],
            cb_pad=CB_PAD,
            cb_size=CB_SIZE,
            cb_label="Equivalent plastic strain / %",
        )
    else:
        axs[1].text(
            0.5,
            0.5,
            r"$\bar{\varepsilon}^{\,p} / \%$ not in VTU",
            ha="center",
            va="center",
            transform=axs[1].transAxes,
        )
        axs[1].set_axis_off()

    vmin, vmax = get_range("strain_tr")
    _plot_contourf_compat(
        meshE,
        fig,
        axs[2],
        strain_tr,
        cmap="RdBu_r",
        vmin=vmin,
        vmax=vmax,
        title=r"$\mathrm{tr}(\boldsymbol{\varepsilon})$",
        fontsize=plt.rcParams["xtick.labelsize"],
        cb_labelsize=plt.rcParams["xtick.labelsize"],
        cb_pad=CB_PAD,
        cb_size=CB_SIZE,
    )

    _savefig(fig, post_dir / "end_fields_tEnd_uRel_epsp_trEps.png")


def compute_global_ranges(
    results: list[Any],
    MASTER: dict,
    *,
    prefix: str | None = None,
    t0: float = 3.0,
    epsp_candidates: tuple[str, ...] = (
        "EquivalentPlasticStrain",
        "eps_p",
        "epsp",
        "equivalent_plastic_strain",
    ),
) -> dict[str, tuple[float, float]]:
    prefix = prefix or MASTER["output"]["prefix"]
    strain_tr = ot.variables.strain.trace

    mins = {"u_rel": np.inf, "strain_tr": np.inf, "epsp": np.inf}
    maxs = {"u_rel": -np.inf, "strain_tr": -np.inf, "epsp": -np.inf}

    for r in results:
        case_dir = Path(_case_attr(r, "case_dir"))
        vtus = list_raw_vtus(case_dir / "results", prefix)
        if not vtus:
            continue

        _, _, vtu_t0 = find_vtu_at_time(vtus, t0)
        vtu_end = vtus[-1][2]

        mesh0 = ot.MeshSeries(str(vtu_t0))[-1]
        meshE = ot.MeshSeries(str(vtu_end))[-1]

        urel = displacement_norm_relative_mm(meshE, mesh0).reshape(-1)
        mins["u_rel"] = float(min(mins["u_rel"], np.nanmin(urel)))
        maxs["u_rel"] = float(max(maxs["u_rel"], np.nanmax(urel)))

        e = np.asarray(strain_tr.transform(meshE)).reshape(-1)
        mins["strain_tr"] = float(min(mins["strain_tr"], np.nanmin(e)))
        maxs["strain_tr"] = float(max(maxs["strain_tr"], np.nanmax(e)))

        pd = meshE.point_data or {}
        cd = meshE.cell_data or {}
        epsp_name = next((k for k in epsp_candidates if (k in pd) or (k in cd)), None)
        if epsp_name:
            ep = _flatten_any(
                pd.get(epsp_name) if epsp_name in pd else cd.get(epsp_name)
            )
            if ep.size:
                mins["epsp"] = float(min(mins["epsp"], np.nanmin(ep)))
                maxs["epsp"] = float(max(maxs["epsp"], np.nanmax(ep)))

    out: dict[str, tuple[float, float]] = {}
    for k, vmin in mins.items():
        vmax = maxs[k]
        if np.isfinite(vmin) and np.isfinite(vmax):
            out[k] = (float(vmin), float(vmax))
    return out


def merge_ranges(
    *,
    global_ranges: dict[str, tuple[float, float]],
    user_ranges: dict[str, tuple[float | None, float | None]] | None,
) -> dict[str, tuple[float | None, float | None]]:
    user_ranges = user_ranges or {}
    keys = ("u_rel", "epsp", "strain_tr")
    merged: dict[str, tuple[float | None, float | None]] = {}
    for key in keys:
        if key in user_ranges:
            merged[key] = user_ranges[key]
        elif key in global_ranges:
            merged[key] = global_ranges[key]
        else:
            merged[key] = (None, None)
    return merged


def run_ssr_postprocessing(
    results: list[Any],
    MASTER: dict,
    *,
    prefix: str | None = None,
    pointA_xyz: np.ndarray | None = None,
    ranges: dict[str, tuple[float | None, float | None]] | None = None,
    auto_compute_missing_ranges: bool = True,
    colorbar_mode: str = "global",
    prefer_ip_stress_for_baseline: bool = True,
) -> None:
    prefix = prefix or MASTER["output"]["prefix"]
    if pointA_xyz is None:
        pointA_xyz = np.array(
            [MASTER["geometry"]["x_crest"], MASTER["geometry"]["H_top"], 0.0],
            dtype=float,
        )

    if colorbar_mode not in ("global", "case"):
        msg = "colorbar_mode must be 'global' or 'case'."
        raise ValueError(msg)

    t0, _t1, _smin, _dt = ssr_params(MASTER)

    user_ranges = ranges or {}
    global_ranges: dict[str, tuple[float, float]] = {}
    if auto_compute_missing_ranges and colorbar_mode == "global":
        global_ranges = compute_global_ranges(results, MASTER, prefix=prefix, t0=t0)

    print("[POST] SSR postprocessing")
    for r in results:
        case_dir = Path(_case_attr(r, "case_dir"))
        name = str(_case_attr(r, "name", case_dir.name))

        results_dir = case_dir / "results"
        post_dir = case_dir / "post"
        post_dir.mkdir(parents=True, exist_ok=True)

        log_path = results_dir / "ogs.log"
        vtus = list_raw_vtus(results_dir, prefix)

        if not log_path.exists():
            print(f"[SKIP] {name}: missing ogs.log at {log_path}")
            continue
        if not vtus:
            print(
                f"[SKIP] {name}: no VTUs found in {results_dir} with prefix '{prefix}'"
            )
            continue

        last_ok_step, t_ok, _, t_fail = parse_last_ok_and_first_fail(log_path)

        # ssr_schedule returns: F, s, c, phi, psi
        F_ok, s_ok, c_ok, phi_ok, psi_ok = ssr_schedule(t_ok, MASTER)
        if t_fail is not None:
            F_fail, *_ = ssr_schedule(t_fail, MASTER)
            fos_text = f"FoS in [{F_ok:.6f}, {F_fail:.6f})"
        else:
            fos_text = f"FoS >= {F_ok:.6f}"

        last_vtu = vtus[-1][2]
        t_end = float(vtus[-1][1])

        _, t0_found, vtu_t0 = find_vtu_at_time(vtus, t0)

        meta = _parse_case_name(name)
        ts_type = _read_time_stepping_type(case_dir / "project.prj")

        print()
        print("=" * 100)
        print("[CASE]")
        if meta is not None:
            print(f"  mesh size h: {meta['h']} m")
            print(f"  element order: {meta['element']} ({meta['ord']})")
            print(f"  minimum dt: {meta['dtmin']}")
            print(f"  nonlinear reltol: {meta['tol']}")
            print(f"  load: {meta['load_tag']} (label: {meta['label']})")
        if ts_type is not None:
            print(f"  time stepping: {ts_type}")
        print(f"  VTU files: {len(vtus)}")
        print(f"  last converged step: {last_ok_step}")
        print(f"  last converged time before failure: {t_ok:.6g}")
        print(
            f"  state at last converged time: s={s_ok:.6g}, c={c_ok:.6g} Pa, "
            f"phi={phi_ok:.6g} deg, psi={psi_ok:.6g} deg"
        )
        print(f"  FoS estimate: {fos_text}")
        print("=" * 100)

        if auto_compute_missing_ranges:
            if colorbar_mode == "global":
                ranges_final = merge_ranges(
                    global_ranges=global_ranges, user_ranges=user_ranges
                )
            else:
                case_ranges = compute_global_ranges([r], MASTER, prefix=prefix, t0=t0)
                ranges_final = merge_ranges(
                    global_ranges=case_ranges, user_ranges=user_ranges
                )
        else:
            ranges_final = merge_ranges(global_ranges={}, user_ranges=user_ranges)

        # build raw time series for point probing + F(t)
        pvd_path = post_dir / f"{name}_raw_series.pvd"
        write_pvd_from_vtus(vtus, pvd_path)
        ms = ot.MeshSeries(pvd_path)

        times = np.array([tt for _, tt, _ in vtus], dtype=float)
        Fvals = np.array([ssr_schedule(t, MASTER)[0] for t in times], dtype=float)

        # time-series plots
        print(
            "[PLOT] time series: F(t) and probe at slope crest "
            f"({pointA_xyz[0]:.6g}, {pointA_xyz[1]:.6g}, {pointA_xyz[2]:.6g})"
        )
        plot_F_vs_time(
            times, Fvals, out_png=post_dir / "F_vs_time.png", t_ok=t_ok, t_fail=t_fail
        )
        plot_point_disp_vs_time(ms, pointA_xyz, out_png=post_dir / "pointA_disp.png")

        # baseline plots (t0)
        print(
            f"[PLOT] baseline stress fields at end of Stage 3 "
            f"(hold/equilibration, before SSR) "
            f"(t={t0_found:.6g}, vtu='{vtu_t0.name}')"
        )
        plot_baseline_stress_row_t0(
            post_dir=post_dir,
            vtu_t0=vtu_t0,
            prefer_ip_stress=prefer_ip_stress_for_baseline,
        )

        # end plots (t_end)
        print(f"[PLOT] end fields @ t_end={t_end:.6g} (vtu='{last_vtu.name}')")
        plot_end_fields_row_tEnd(
            post_dir=post_dir,
            vtu_t0=vtu_t0,
            vtu_end=last_vtu,
            ranges=ranges_final,
        )

    print("[POST] done.")
