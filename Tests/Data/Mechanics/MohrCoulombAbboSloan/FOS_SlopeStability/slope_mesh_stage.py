from __future__ import annotations

import math
import shutil
import traceback
import xml.etree.ElementTree as ET
from collections.abc import Sequence
from dataclasses import dataclass
from itertools import pairwise
from pathlib import Path
from typing import Literal

import gmsh
import matplotlib.pyplot as plt
import meshio
import numpy as np
import ogstools as ot
from slope_utils import case_folder_name

# ============================================================
# PRJ inspection helpers (we NEVER change called mesh names)
# ============================================================


def _all_vtu_mesh_names_in_prj(template_prj: Path) -> set[str]:
    """
    Collect all *.vtu filenames referenced inside ANY <mesh> tag.
    This includes bulk mesh + boundary meshes used in BCs, submesh output, etc.
    """
    root = ET.parse(template_prj).getroot()
    names: set[str] = set()

    # Any <mesh> tag anywhere
    for node in root.findall(".//mesh"):
        if node.text:
            t = node.text.strip()
            if t.lower().endswith(".vtu"):
                names.add(t)

    # <meshes><mesh>...</mesh></meshes> (redundant but safe)
    for node in root.findall(".//meshes/mesh"):
        if node.text:
            t = node.text.strip()
            if t.lower().endswith(".vtu"):
                names.add(t)

    return names


def read_template_bulk_mesh_filename(template_prj: Path) -> str:
    root = ET.parse(template_prj).getroot()
    node = root.find(".//mesh")
    if node is None or not node.text or not node.text.strip():
        msg = f"Template PRJ has no <mesh>...</mesh>: {template_prj}"
        raise RuntimeError(msg)
    return node.text.strip()


# ============================================================
# Mesh generation (Gmsh) -> .msh
# ============================================================

LoadMode = Literal["all", "from_slope_crest", "from_right_edge", "top_segments"]


def _dist(a: tuple[float, float], b: tuple[float, float]) -> float:
    return math.hypot(b[0] - a[0], b[1] - a[1])


def _polyline_total_length(xy: list[tuple[float, float]]) -> float:
    return sum(_dist(xy[i], xy[i + 1]) for i in range(len(xy) - 1))


def _interp(
    a: tuple[float, float], b: tuple[float, float], t: float
) -> tuple[float, float]:
    return (a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1]))


def _insert_points_by_s(
    xy: list[tuple[float, float]],
    s_cuts: list[float],
    tol: float = 1e-12,
) -> list[tuple[float, float]]:
    L = _polyline_total_length(xy)
    if L <= 0:
        msg = "Polyline has zero length."
        raise ValueError(msg)

    s_clean = sorted({s for s in s_cuts if (tol < s < L - tol)})
    if not s_clean:
        return xy[:]

    new_xy: list[tuple[float, float]] = [xy[0]]
    s_idx = 0
    s_target = s_clean[s_idx]
    s_acc = 0.0

    for i in range(len(xy) - 1):
        a, b = xy[i], xy[i + 1]
        seg = _dist(a, b)
        if seg <= tol:
            continue

        while s_idx < len(s_clean) and s_acc + seg >= s_target - tol:
            t = (s_target - s_acc) / seg
            t = min(1.0, max(0.0, t))
            pcut = _interp(a, b, t)
            if _dist(new_xy[-1], pcut) > tol:
                new_xy.append(pcut)

            s_idx += 1
            if s_idx >= len(s_clean):
                break
            s_target = s_clean[s_idx]

        if _dist(new_xy[-1], b) > tol:
            new_xy.append(b)

        s_acc += seg

    return new_xy


def _select_top_lines_overlapping(
    top_lines: list[int],
    top_xy: list[tuple[float, float]],
    a: float,
    b: float,
    tol: float = 1e-12,
) -> list[int]:
    s_breaks = [0.0]
    for i in range(len(top_xy) - 1):
        s_breaks.append(s_breaks[-1] + _dist(top_xy[i], top_xy[i + 1]))

    out: list[int] = []
    for k, lt in enumerate(top_lines):
        s0, s1 = s_breaks[k], s_breaks[k + 1]
        if (s1 > a + tol) and (s0 < b - tol):
            out.append(lt)
    return out


def create_slope_msh(msh_path: Path, GEOMETRY: dict, LOAD: dict, MESH: dict) -> Path:
    msh_path = Path(msh_path)
    model_name = msh_path.stem

    Lb = float(GEOMETRY["L_bottom"])
    Ht = float(GEOMETRY["H_top"])
    Hb = float(GEOMETRY["H_bench"])
    x_toe = float(GEOMETRY["x_toe"])
    top_breaks = [float(x) for x in GEOMETRY["top_breaks"]]
    h = float(GEOMETRY["h"])

    p0_xy = (0.0, 0.0)
    p1_xy = (Lb, 0.0)
    top_xy_base = [(x, Ht) for x in top_breaks]
    p6_xy = (x_toe, Hb)
    p7_xy = (0.0, Hb)

    mode: LoadMode = str(LOAD.get("mode", "all")).strip().lower()  # type: ignore[assignment]
    clamp = bool(LOAD.get("clamp", True))

    L_top = _polyline_total_length(top_xy_base)
    s_cuts: list[float] = []

    if mode in ("from_slope_crest", "from_right_edge"):
        L_load = float(LOAD["L_load"])
        L_load_eff = min(L_load, L_top) if clamp else L_load
        s_cuts = [L_load_eff] if mode == "from_right_edge" else [L_top - L_load_eff]
    elif mode == "top_segments":
        L1 = float(LOAD["L_top_1"])
        L2 = float(LOAD["L_top_2"])
        s_cuts = [L1, L1 + L2]

    top_xy = _insert_points_by_s(top_xy_base, s_cuts)

    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Verbosity", 1)
        gmsh.model.add(model_name)
        geo = gmsh.model.geo

        tag = 1
        p0 = geo.addPoint(*p0_xy, 0.0, h, tag=tag)
        tag += 1
        p1 = geo.addPoint(*p1_xy, 0.0, h, tag=tag)
        tag += 1

        top_pt_tags: list[int] = []
        for xy in top_xy:
            pt = geo.addPoint(xy[0], xy[1], 0.0, h, tag=tag)
            top_pt_tags.append(pt)
            tag += 1

        p6 = geo.addPoint(*p6_xy, 0.0, h, tag=tag)
        tag += 1
        p7 = geo.addPoint(*p7_xy, 0.0, h, tag=tag)
        tag += 1

        l_bottom = geo.addLine(p0, p1)
        l_right = geo.addLine(p1, top_pt_tags[0])

        top_lines: list[int] = []
        for a_pt, b_pt in pairwise(top_pt_tags):
            top_lines.append(geo.addLine(a_pt, b_pt))

        l_slope = geo.addLine(top_pt_tags[-1], p6)
        l_bench = geo.addLine(p6, p7)
        l_left = geo.addLine(p7, p0)

        loop_edges = [l_bottom, l_right] + top_lines + [l_slope, l_bench, l_left]
        cl = geo.addCurveLoop(loop_edges)
        srf = geo.addPlaneSurface([cl])

        geo.synchronize()

        gmsh.option.setNumber("Mesh.Algorithm", int(MESH.get("algo_2d", 6)))
        gmsh.option.setNumber("Mesh.ElementOrder", 1)
        if bool(MESH.get("recombine", False)):
            gmsh.option.setNumber("Mesh.RecombineAll", 1)

        pg = gmsh.model.addPhysicalGroup
        pn = gmsh.model.setPhysicalName

        pn(1, pg(1, [l_bottom]), "bottom")
        pn(1, pg(1, [l_right]), "right")
        pn(1, pg(1, [l_left]), "left")
        pn(1, pg(1, [l_slope]), "slope")
        pn(1, pg(1, [l_bench]), "bench")
        pn(1, pg(1, top_lines), "top_all")

        top_load: list[int] = []
        top_rest: list[int] = []

        if mode == "all":
            top_load = top_lines[:]

        elif mode == "from_right_edge":
            L_load_eff = (
                min(float(LOAD["L_load"]), L_top) if clamp else float(LOAD["L_load"])
            )
            top_load = _select_top_lines_overlapping(top_lines, top_xy, 0.0, L_load_eff)
            top_load_set = set(top_load)
            top_rest = [lt for lt in top_lines if lt not in top_load_set]

        elif mode == "from_slope_crest":
            L_load_eff = (
                min(float(LOAD["L_load"]), L_top) if clamp else float(LOAD["L_load"])
            )
            a, b = L_top - L_load_eff, L_top
            top_load = _select_top_lines_overlapping(top_lines, top_xy, a, b)
            top_load_set = set(top_load)
            top_rest = [lt for lt in top_lines if lt not in top_load_set]

        elif mode == "top_segments":
            L1 = float(LOAD["L_top_1"])
            L2 = float(LOAD["L_top_2"])
            s1, s2 = L1, L1 + L2

            z1 = _select_top_lines_overlapping(top_lines, top_xy, 0.0, s1)
            z2 = _select_top_lines_overlapping(top_lines, top_xy, s1, s2)
            z3 = _select_top_lines_overlapping(top_lines, top_xy, s2, L_top)

            pn(1, pg(1, z1), "top_zone_1")
            pn(1, pg(1, z2), "top_zone_2")
            pn(1, pg(1, z3), "top_zone_3")

            zones = {int(z) for z in LOAD.get("loaded_zones", [2])}
            top_load_set: set[int] = set()
            if 1 in zones:
                top_load_set.update(z1)
            if 2 in zones:
                top_load_set.update(z2)
            if 3 in zones:
                top_load_set.update(z3)

            top_load = [lt for lt in top_lines if lt in top_load_set]
            top_rest = [lt for lt in top_lines if lt not in top_load_set]

        pn(1, pg(1, top_load), "top_load")
        if top_rest:
            pn(1, pg(1, top_rest), "top_rest")

        pn(2, pg(2, [srf]), "domain")

        if bool(MESH.get("generate_mesh", True)):
            gmsh.model.mesh.generate(2)

        msh_path.parent.mkdir(parents=True, exist_ok=True)
        gmsh.write(str(msh_path))
        return msh_path
    finally:
        gmsh.finalize()


# ============================================================
# MSH -> VTU (NO renaming; bulk uses template PRJ filename)
# ============================================================


def msh_to_vtus_for_template_prj(
    msh_file: Path,
    mesh_dir: Path,
    template_prj: Path,
    *,
    overwrite: bool = True,
) -> dict[str, Path]:
    mesh_dir = Path(mesh_dir)
    mesh_dir.mkdir(parents=True, exist_ok=True)

    meshes = ot.Meshes.from_gmsh(msh_file, log=False)
    if "domain" not in meshes:
        msg = "No domain mesh found in .msh"
        raise RuntimeError(msg)

    bulk_mesh = meshes["domain"]
    bulk_filename = read_template_bulk_mesh_filename(template_prj)

    written: dict[str, Path] = {}

    bulk_path = mesh_dir / bulk_filename
    if overwrite or not bulk_path.exists():
        ot.mesh.save(bulk_mesh, str(bulk_path))
    written["__bulk__"] = bulk_path

    for name, m in meshes.items():
        if name == "domain":
            continue
        p = mesh_dir / f"{name}.vtu"
        if overwrite or not p.exists():
            ot.mesh.save(m, str(p))
        written[name] = p

    expected = _all_vtu_mesh_names_in_prj(template_prj)
    missing = [fn for fn in sorted(expected) if not (mesh_dir / fn).exists()]
    if missing:
        msg = (
            "Template PRJ expects meshes that are not present in this case mesh folder:\n"
            + "\n".join(f"  - {fn}" for fn in missing)
        )
        raise FileNotFoundError(msg)

    return written


# ============================================================
# Quadratic conversion (in-place, filenames unchanged)
# ============================================================


def identify_subdomains_in_place(
    mesh_dir: Path, template_prj: Path, *, eps: float = 1e-14
) -> None:
    mesh_dir = Path(mesh_dir)

    bulk_fn = read_template_bulk_mesh_filename(template_prj)
    bulk_path = mesh_dir / bulk_fn

    expected = _all_vtu_mesh_names_in_prj(template_prj)
    boundary_files = [mesh_dir / fn for fn in sorted(expected) if fn != bulk_fn]

    cli = ot.cli()
    cli.identifySubdomains(
        "-m",
        str(bulk_path),
        f"-o {mesh_dir}/",
        "-f",
        "-s",
        str(eps),
        "--",
        *[str(p) for p in boundary_files],
    )


def make_quadratic_in_place(mesh_dir: Path, template_prj: Path) -> None:
    cli = ot.cli()
    bulk_vtu = Path(mesh_dir) / read_template_bulk_mesh_filename(template_prj)
    for p in sorted(Path(mesh_dir).glob("*.vtu")):
        if p == bulk_vtu:
            # Ensure the process mesh is strictly 2D before/after conversion.
            keep_2d_cells_only(p)

        m = meshio.read(str(p))
        cell_types = {cb.type for cb in m.cells}
        if cell_types & {"triangle6", "quad8", "quad9", "line3"}:
            continue
        if cell_types & {"triangle", "quad", "line"}:
            cli.createQuadraticMesh("-i", str(p), "-o", str(p))

        if p == bulk_vtu:
            keep_2d_cells_only(p)


_TWO_D_CELL_TYPES = ("triangle", "triangle6", "quad", "quad8", "quad9")


def _filtered_cells_and_cell_data(
    vtu: meshio.Mesh, allowed_cell_types: tuple[str, ...]
) -> tuple[list[meshio.CellBlock], dict[str, list[np.ndarray]]]:
    # Keep cells and cell_data block indices aligned after type-based filtering.
    new_cells: list[meshio.CellBlock] = []
    new_cell_data: dict[str, list[np.ndarray]] = {}
    for cb_i, cb in enumerate(vtu.cells):
        if cb.type not in allowed_cell_types:
            continue
        new_cells.append(meshio.CellBlock(cb.type, cb.data))
        for name, blocks in vtu.cell_data.items():
            if cb_i < len(blocks):
                new_cell_data.setdefault(name, []).append(blocks[cb_i])
    return new_cells, new_cell_data


def strip_material_ids_in_place(mesh_dir: Path) -> None:
    for p in sorted(Path(mesh_dir).glob("*.vtu")):
        m = meshio.read(str(p))
        if m.cell_data:
            m.cell_data = {
                k: v for k, v in m.cell_data.items() if "material" not in k.lower()
            }
        if m.point_data:
            m.point_data = {
                k: v for k, v in m.point_data.items() if "material" not in k.lower()
            }
        if m.field_data:
            m.field_data = {
                k: v for k, v in m.field_data.items() if "material" not in k.lower()
            }
        meshio.write(str(p), m)


def assign_material_ids_all_zero(bulk_vtu: Path) -> None:
    vtu = meshio.read(str(bulk_vtu))
    new_cells, _ = _filtered_cells_and_cell_data(vtu, _TWO_D_CELL_TYPES)
    new_cell_data: dict[str, list[np.ndarray]] = {"MaterialIDs": []}

    for cb in new_cells:
        mats = np.zeros(len(cb.data), dtype=np.int32)
        new_cell_data["MaterialIDs"].append(mats)

    if not new_cells:
        msg = "No 2D cells found to assign MaterialIDs."
        raise ValueError(msg)

    out = meshio.Mesh(points=vtu.points, cells=new_cells, cell_data=new_cell_data)
    meshio.write(str(bulk_vtu), out)


def keep_2d_cells_only(bulk_vtu: Path) -> None:
    vtu = meshio.read(str(bulk_vtu))
    new_cells, new_cell_data = _filtered_cells_and_cell_data(vtu, _TWO_D_CELL_TYPES)
    if not new_cells:
        msg = "Bulk mesh has no 2D cells after filtering."
        raise ValueError(msg)
    out = meshio.Mesh(points=vtu.points, cells=new_cells, cell_data=new_cell_data)
    meshio.write(str(bulk_vtu), out)


# ============================================================
# Mesh plotting (matplotlib only; saves PNG)
# ============================================================


def plot_mesh_from_msh(
    msh_path: Path,
    out_png: Path,
    *,
    show: bool = True,
    close: bool = False,
    domain_alpha: float = 0.20,
    boundary_lw: float = 3.0,
    show_boundaries: tuple[str, ...] = (
        "slope",
        "bottom",
        "top_load",
        "top_rest",
        "left",
        "right",
        "bench",
    ),
    arrow_count: int = 7,
    debug: bool = False,
):
    msh_path = Path(msh_path)
    out_png = Path(out_png)
    out_png.parent.mkdir(parents=True, exist_ok=True)

    m = meshio.read(msh_path)
    pts = np.asarray(m.points)[:, :2]

    def to_tris(celltype: str, conn: np.ndarray) -> np.ndarray:
        if celltype == "triangle":
            return conn
        if celltype == "triangle6":
            return conn[:, :3]
        if celltype in ("quad", "quad8", "quad9"):
            q = conn[:, :4]
            return np.vstack([q[:, [0, 1, 2]], q[:, [0, 2, 3]]])
        return np.zeros((0, 3), dtype=int)

    tris_list: list[np.ndarray] = []
    for cb in m.cells:
        if cb.type not in ("triangle", "triangle6", "quad", "quad8", "quad9"):
            continue
        conn = np.asarray(cb.data, dtype=int)
        t = to_tris(cb.type, conn)
        if t.size:
            tris_list.append(t)

    if not tris_list:
        msg = "No 2D cells found in msh for plotting."
        raise RuntimeError(msg)
    tris = np.vstack(tris_list)

    fig, ax = plt.subplots(figsize=(10, 3.2), dpi=150)
    ax.triplot(pts[:, 0], pts[:, 1], tris, linewidth=0.25, alpha=domain_alpha)

    meshes1 = ot.Meshes.from_gmsh(msh_path, dim=[1], log=False)

    def _pts2(mesh_obj) -> np.ndarray:
        p = mesh_obj.points
        return np.asarray(p)[:, :2]

    handles, labels = [], []
    cmap = plt.get_cmap("tab10")
    ci = 0
    top_load_pts: np.ndarray | None = None

    for bname in show_boundaries:
        key = next((k for k in meshes1 if bname.lower() in k.lower()), None)
        if key is None:
            continue

        bp = _pts2(meshes1[key])
        if bp.shape[0] < 2:
            continue

        idx = (
            np.argsort(bp[:, 0])
            if np.ptp(bp[:, 0]) >= np.ptp(bp[:, 1])
            else np.argsort(bp[:, 1])
        )

        color = cmap(ci % 10)
        ci += 1
        (ln,) = ax.plot(bp[idx, 0], bp[idx, 1], lw=boundary_lw, alpha=0.95, color=color)
        handles.append(ln)
        labels.append(bname)

        if bname.lower() == "top_load":
            top_load_pts = bp[idx]

    if top_load_pts is not None:
        x = top_load_pts[:, 0]
        y = top_load_pts[:, 1]
        xuniq, uidx = np.unique(x, return_index=True)
        yuniq = y[uidx]

        if xuniq.size >= 2 and arrow_count > 0 and (xuniq.max() - xuniq.min()) > 0:
            xs = np.linspace(float(xuniq.min()), float(xuniq.max()), arrow_count)
            ys = np.interp(xs, xuniq, yuniq)

            y0, y1 = ax.get_ylim()
            arrow_len = 1.5
            ax.set_ylim(y0, y1 + 1.5 * arrow_len)

            for xi, yi in zip(xs, ys, strict=False):
                ax.annotate(
                    "",
                    xy=(xi, yi),
                    xytext=(xi, yi + arrow_len),
                    arrowprops={"arrowstyle": "-|>", "lw": 1.6, "color": "k"},
                    annotation_clip=False,
                )
    elif debug:
        print("No top_load boundary found. Available 1D meshes:", list(meshes1.keys()))

    ax.set_aspect("equal")
    ax.set_xlabel("x / m", fontsize=11)
    ax.set_ylabel("y / m", fontsize=11)
    ax.tick_params(labelsize=10)

    if handles:
        fig.subplots_adjust(right=0.78)
        ax.legend(
            handles,
            labels,
            fontsize=8,
            loc="center left",
            bbox_to_anchor=(1.02, 0.5),
            frameon=True,
        )

    fig.tight_layout()
    fig.savefig(out_png, bbox_inches="tight")

    if show:
        plt.show()
    if close:
        plt.close(fig)

    return fig, ax


# ============================================================
# MeshStage class
# ============================================================

MeshMode = Literal["generate", "copy", "use"]


def _check_mesh_folder_against_prj(mesh_dir: Path, template_prj: Path) -> None:
    """
    Ensures *all* VTUs referenced by the PRJ exist in mesh_dir.
    This is the key sanity check for copy/use paths.
    """
    mesh_dir = Path(mesh_dir)
    if not mesh_dir.is_dir():
        msg = f"mesh_dir is not a directory: {mesh_dir}"
        raise FileNotFoundError(msg)

    expected = _all_vtu_mesh_names_in_prj(template_prj)
    missing = [fn for fn in sorted(expected) if not (mesh_dir / fn).exists()]
    if missing:
        msg = (
            "Mesh folder is missing VTU files required by template PRJ:\n"
            + "\n".join(f"  - {fn}" for fn in missing)
            + f"\n\nmesh_dir = {mesh_dir.resolve()}"
        )
        raise FileNotFoundError(msg)


def _has_all_vtus(mesh_dir: Path, template_prj: Path) -> bool:
    expected = _all_vtu_mesh_names_in_prj(template_prj)
    return all((Path(mesh_dir) / fn).exists() for fn in expected)


def _copy_tree_clean(src: Path, dst: Path, *, overwrite: bool) -> None:
    """
    Copy src -> dst. If overwrite=True, delete dst first for a clean copy.
    """
    src = Path(src)
    dst = Path(dst)
    if not src.is_dir():
        msg = f"Source mesh dir not found: {src}"
        raise FileNotFoundError(msg)

    if overwrite and dst.exists():
        shutil.rmtree(dst)
    dst.mkdir(parents=True, exist_ok=True)

    for p in src.rglob("*"):
        rel = p.relative_to(src)
        out = dst / rel
        if p.is_dir():
            out.mkdir(parents=True, exist_ok=True)
        else:
            out.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(p, out)


@dataclass
class CaseSpec:
    label: str
    h: float
    element_order: int
    nonlinear_reltol: float
    minimum_dt: float
    load: dict
    load_tag: str

    mesh_mode: MeshMode = "generate"
    src_mesh_dir: Path | None = None  # required for copy/use

    regenerate: bool = True  # only used when mesh_mode="generate"
    copy_overwrite: bool = True  # only used when mesh_mode="copy"
    plot_mesh: bool = True
    mesh_dir: Path | None = None


@dataclass
class MeshStageResult:
    name: str
    case_dir: Path
    ok: bool
    error: str | None = None
    msh: Path | None = None
    mesh_dir: Path | None = None
    mesh_mode: MeshMode | None = None
    mesh_source: Path | None = None


class MeshStage:
    def __init__(self, *, out_root: Path) -> None:
        self.out_root = Path(out_root)
        self.out_root.mkdir(parents=True, exist_ok=True)

    def build_one(
        self, *, basename: str, template_prj: Path, master: dict, case: CaseSpec
    ) -> MeshStageResult:
        case_name = case_folder_name(
            basename,
            element_order=case.element_order,
            h=case.h,
            dtmin=case.minimum_dt,
            tol=case.nonlinear_reltol,
            load_tag=case.load_tag,
            label=case.label,
        )

        case_dir = self.out_root / case_name
        case_mesh_dir = case_dir / "mesh"
        results_dir = case_dir / "results"
        case_mesh_dir.mkdir(parents=True, exist_ok=True)
        results_dir.mkdir(parents=True, exist_ok=True)

        msh_path = case_mesh_dir / f"{master['mesh_name']}.msh"

        try:
            print(f"\n[MESH-CASE] {case_name}")
            print(f"  mesh_mode = {case.mesh_mode}")

            # -----------------------------------------
            # Mode: USE existing mesh folder in-place
            # -----------------------------------------
            if case.mesh_mode == "use":
                if case.src_mesh_dir is None:
                    msg = "mesh_mode='use' requires src_mesh_dir"
                    raise ValueError(msg)
                src = Path(case.src_mesh_dir)
                _check_mesh_folder_against_prj(src, template_prj)
                return MeshStageResult(
                    name=case_name,
                    case_dir=case_dir,
                    ok=True,
                    msh=None,
                    mesh_dir=src,
                    mesh_mode=case.mesh_mode,
                    mesh_source=src,
                )

            # -----------------------------------------
            # Mode: COPY existing mesh folder into case
            # -----------------------------------------
            if case.mesh_mode == "copy":
                if case.src_mesh_dir is None:
                    msg = "mesh_mode='copy' requires src_mesh_dir"
                    raise ValueError(msg)
                src = Path(case.src_mesh_dir)
                _check_mesh_folder_against_prj(src, template_prj)
                _copy_tree_clean(src, case_mesh_dir, overwrite=case.copy_overwrite)
                _check_mesh_folder_against_prj(case_mesh_dir, template_prj)

                return MeshStageResult(
                    name=case_name,
                    case_dir=case_dir,
                    ok=True,
                    msh=None,
                    mesh_dir=case_mesh_dir,
                    mesh_mode=case.mesh_mode,
                    mesh_source=src,
                )

            # -----------------------------------------
            # Mode: GENERATE (gmsh -> vtu)
            # -----------------------------------------
            if case.mesh_mode != "generate":
                msg = f"Unknown mesh_mode: {case.mesh_mode}"
                raise ValueError(msg)

            # If regenerate is False and we already have ALL required VTUs, reuse them
            if (not case.regenerate) and _has_all_vtus(case_mesh_dir, template_prj):
                print("  [MESH] reuse existing case mesh folder (regenerate=False)")
                return MeshStageResult(
                    name=case_name,
                    case_dir=case_dir,
                    ok=True,
                    msh=msh_path if msh_path.exists() else None,
                    mesh_dir=case_mesh_dir,
                    mesh_mode=case.mesh_mode,
                    mesh_source=case_mesh_dir,
                )

            geom = dict(master["geometry"])
            geom["h"] = float(case.h)

            create_slope_msh(
                msh_path, GEOMETRY=geom, LOAD=case.load, MESH=master["mesh"]
            )
            msh_to_vtus_for_template_prj(
                msh_path, case_mesh_dir, template_prj, overwrite=True
            )

            if case.element_order == 2:
                strip_material_ids_in_place(case_mesh_dir)
                make_quadratic_in_place(case_mesh_dir, template_prj)

            identify_subdomains_in_place(case_mesh_dir, template_prj)
            bulk_vtu = case_mesh_dir / read_template_bulk_mesh_filename(template_prj)
            assign_material_ids_all_zero(bulk_vtu)
            keep_2d_cells_only(bulk_vtu)
            _check_mesh_folder_against_prj(case_mesh_dir, template_prj)

            if case.plot_mesh:
                png = case_mesh_dir / "mesh.png"
                plot_mesh_from_msh(msh_path, png, show=True, close=False)
                print(f"[PLOT] saved -> {png.resolve()}")

            return MeshStageResult(
                name=case_name,
                case_dir=case_dir,
                ok=True,
                msh=msh_path,
                mesh_dir=case_mesh_dir,
                mesh_mode=case.mesh_mode,
                mesh_source=case_mesh_dir,
            )

        except Exception as e:
            (case_dir / "mesh_error_trace.txt").write_text(traceback.format_exc())
            return MeshStageResult(
                name=case_name,
                case_dir=case_dir,
                ok=False,
                error=f"{type(e).__name__}: {e}",
                msh=msh_path if msh_path.exists() else None,
                mesh_dir=case_mesh_dir if case_mesh_dir.exists() else None,
                mesh_mode=case.mesh_mode,
                mesh_source=Path(case.src_mesh_dir) if case.src_mesh_dir else None,
            )

    def build_all(
        self,
        *,
        basename: str,
        template_prj: Path,
        master: dict,
        cases: Sequence[CaseSpec],
    ) -> list[MeshStageResult]:
        out: list[MeshStageResult] = []
        for c in cases:
            r = self.build_one(
                basename=basename, template_prj=template_prj, master=master, case=c
            )
            print("[OK]" if r.ok else "[FAILED]", r.name)
            if not r.ok:
                print("  error:", r.error)
            out.append(r)
        return out
