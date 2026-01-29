from __future__ import annotations

import traceback
import xml.etree.ElementTree as ET
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import ogstools as ot
from slope_utils import (
    case_folder_name,
    rho_sr_expr,
    ssr_params,
    ssr_strength_expr,
    stage_params,
    top_load_expr,
)

# ============================================================
# PRJ inspection helpers (for sanity checks)
# ============================================================


def _all_vtu_mesh_names_in_prj(template_prj: Path) -> set[str]:
    root = ET.parse(template_prj).getroot()
    names: set[str] = set()

    for node in root.findall(".//mesh"):
        if node.text:
            t = node.text.strip()
            if t.lower().endswith(".vtu"):
                names.add(t)

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
# Data models
# ============================================================


@dataclass
class CaseSpec:
    label: str
    h: float
    element_order: int
    nonlinear_reltol: float
    minimum_dt: float
    load: dict
    load_tag: str
    plot_mesh: bool = True  # unused by runner, kept for compatibility
    fixed_output_times: tuple[float, ...] = field(default_factory=tuple)
    mesh_dir: Path | None = None


@dataclass
class OGSResult:
    name: str
    case_dir: Path
    ok: bool
    stage: str
    error: str | None = None
    logfile: Path | None = None
    log_tail: str | None = None


# ============================================================
# PRJ update
# ============================================================


def apply_cfg_to_prj(
    *,
    template_prj: Path,
    prj_out: Path,
    master: dict,
    case: CaseSpec,
) -> Path:
    prj = ot.Project(input_file=template_prj, output_file=prj_out)

    def _fmt(val):
        if isinstance(val, bool):
            return str(val)
        return f"{val:.16g}" if isinstance(val, int | float) else str(val)

    def _set(xpath: str, val, *, required: bool = True):
        n = prj.replace_text(_fmt(val), xpath=xpath)
        if n == 0 and required:
            msg = f"No nodes matched xpath: {xpath}"
            raise KeyError(msg)

    def _set_value(param_name: str, val):
        if isinstance(val, str):
            expr = val.strip()
            prj.replace_text(
                "Function",
                xpath=f".//parameters/parameter[name='{param_name}']/type",
            )
            n = prj.replace_text(
                expr,
                xpath=f".//parameters/parameter[name='{param_name}']/expression",
            )
            if n == 0:
                msg = f"Parameter '{param_name}' has no <expression> in template PRJ."
                raise KeyError(msg)
            return

        # numeric -> Constant + <value>/<values>
        sval = _fmt(val)
        prj.replace_text(
            "Constant",
            xpath=f".//parameters/parameter[name='{param_name}']/type",
        )
        n = prj.replace_text(
            sval, xpath=f".//parameters/parameter[name='{param_name}']/value"
        )
        if n == 0:
            n = prj.replace_text(
                sval, xpath=f".//parameters/parameter[name='{param_name}']/values"
            )
        if n == 0:
            msg = f"Parameter '{param_name}' has no <value> or <values> node in template PRJ."
            raise KeyError(msg)

    def _set_values(param_name: str, vec):
        txt = " ".join(_fmt(float(v)) for v in vec)
        n = prj.replace_text(
            txt, xpath=f".//parameters/parameter[name='{param_name}']/values"
        )
        if n == 0:
            n = prj.replace_text(
                txt, xpath=f".//parameters/parameter[name='{param_name}']/value"
            )
        if n == 0:
            msg = f"Parameter '{param_name}' has no <value> or <values> node in template PRJ."
            raise KeyError(msg)

    # ---------------------------------------------------------
    # fixed output times: SSR start time + one extra time (case-specific)
    t0 = float(master.get("ssr", {}).get("t0", 3.0))
    times = [t0]
    times.append(4.7)
    times_txt = " ".join(_fmt(x) for x in times)
    _set(".//time_loop/output/fixed_output_times", times_txt)
    # ---------------------------------------------------------
    # Build effective scalar params WITHOUT mutating master
    # ---------------------------------------------------------
    scalars = master["material"]["parameters"]["scalar"]
    scalars_eff = dict(scalars)  # copy

    # ---- staging expressions (rho_sr, top_load)
    rho0 = float(scalars_eff["rho_sr0"])
    top0 = float(scalars_eff["top_load0"])
    scalars_eff["rho_sr"] = rho_sr_expr(master, rho0=rho0)
    scalars_eff["top_load"] = top_load_expr(master, top0=top0)

    # ---- SSR expressions (Cohesion, FrictionAngle, DilatancyAngle)
    ssr = master.get("ssr", {})
    if ssr:
        s_expr = ssr_strength_expr(master)

        cohesion0 = float(scalars_eff["cohesion_rs0"])
        phi0 = float(scalars_eff["FrictionAngle_rs0"])
        psi0 = float(scalars_eff["DilatancyAngle_rs0"])

        scalars_eff["Cohesion"] = f"{cohesion0:.16g} * {s_expr}"
        scalars_eff["FrictionAngle"] = (
            f"atan(({s_expr}) * tan({phi0:.16g}*pi/180.0)) * 180.0/pi"
        )
        scalars_eff["DilatancyAngle"] = (
            f"atan(({s_expr}) * tan({psi0:.16g}*pi/180.0)) * 180.0/pi"
        )

    # ---------------------------------------------------------
    # Write parameters
    # ---------------------------------------------------------
    for name, val in scalars_eff.items():
        _set_value(name, val)

    for name, vec in master["material"]["parameters"]["vector"].items():
        _set_values(name, vec)

    # time stepping
    ts = dict(master["time_stepping"])
    ts["minimum_dt"] = float(case.minimum_dt)
    for node, val in ts.items():
        _set(f".//time_loop/processes/process/time_stepping/{node}", val)

    # solver reltol
    _set(
        ".//time_loop/processes/process/convergence_criterion/reltol",
        float(case.nonlinear_reltol),
    )

    # element orders
    displacement_order = 1 if case.element_order == 1 else 2
    integration_order = 2 if case.element_order == 1 else 4
    _set(
        ".//process_variables/process_variable[name='displacement']/order",
        displacement_order,
    )
    _set(".//processes/process/integration_order", integration_order)

    # output
    _set(".//time_loop/output/prefix", master["output"]["prefix"])
    _set(".//time_loop/output/suffix", master["output"]["suffix"])

    # body force
    if "physics" in master and "specific_body_force" in master["physics"]:
        bf = master["physics"]["specific_body_force"]
        _set(
            ".//processes/process/specific_body_force",
            " ".join(_fmt(float(x)) for x in bf),
        )

    prj.write_input()
    return prj_out


# ============================================================
# Plotting
# ============================================================
def plot_staging_and_ssr_curves(
    *,
    master: dict,
    case_dir: Path,
    _case_name: str,  # unused by design (kept for compatibility), satisfies ARG001
    t_end: float | None = None,
    n: int = 2001,
    show: bool = False,
) -> tuple[Path, Path]:
    """
    Plot the *actual* staging + SSR schedules implied by MASTER (function-based).
    Saves PNG + PDF into case_dir/results.
    If show=True, also displays the figure before running OGS.
    """
    case_dir = Path(case_dir)
    out_dir = case_dir / "results"
    out_dir.mkdir(parents=True, exist_ok=True)

    scalars = master["material"]["parameters"]["scalar"]

    # Time horizon
    if t_end is None:
        t_end = float(master["time_stepping"]["t_end"])
    t = np.linspace(0.0, float(t_end), int(n))

    # Baselines
    rho0 = float(scalars["rho_sr0"])
    q0 = float(scalars["top_load0"])
    c0 = float(scalars["cohesion_rs0"])
    phi0_deg = float(scalars["FrictionAngle_rs0"])
    psi0_deg = float(scalars["DilatancyAngle_rs0"])

    # Stage times
    rho_t1, load_t0, load_t1 = stage_params(master)

    # SSR times
    t0, t1, smin, dt = ssr_params(master)

    # -------------------------
    # Schedules (function-based)
    # -------------------------
    rho_ratio = np.clip(t / rho_t1, 0.0, 1.0)
    q_ratio = np.clip((t - load_t0) / (load_t1 - load_t0), 0.0, 1.0)

    s_lin = 1.0 - (1.0 - smin) * (t - t0) / dt
    s_t = np.clip(s_lin, smin, 1.0)
    s_t = np.where(t < t0, 1.0, s_t)

    c_ratio = s_t

    phi0 = np.deg2rad(phi0_deg)
    psi0 = np.deg2rad(psi0_deg)
    phi_t = np.rad2deg(np.arctan(np.tan(phi0) * s_t))
    psi_t = np.rad2deg(np.arctan(np.tan(psi0) * s_t))
    phi_ratio = phi_t / phi0_deg if phi0_deg != 0 else np.nan * np.ones_like(t)
    psi_ratio = psi_t / psi0_deg if psi0_deg != 0 else np.nan * np.ones_like(t)

    F_t = 1.0 / np.clip(s_t, 1e-12, None)
    t_F_end = 8.0
    F_plot = F_t.copy()
    F_plot[(t < t0) | (t > t_F_end)] = np.nan

    # -------------------------
    # Plot styling
    # -------------------------
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["DejaVu Serif"],
            "font.size": 16,
            "axes.labelsize": 16,
            "axes.titlesize": 16,
            "legend.fontsize": 12,
            "xtick.labelsize": 14,
            "ytick.labelsize": 14,
            "lines.linewidth": 2.2,
            "axes.linewidth": 1.0,
            "grid.linewidth": 0.7,
            "savefig.dpi": 300,
            "mathtext.fontset": "stix",
        }
    )

    fig, ax = plt.subplots(figsize=(13.2, 5.8))

    # Main normalized curves
    (l_rho,) = ax.plot(t, rho_ratio, label=r"$\dfrac{\rho(t)}{\rho_0}$")
    (l_q,) = ax.plot(t, q_ratio, label=r"$\dfrac{q(t)}{q_0}$")
    (l_c,) = ax.plot(t, c_ratio, label=r"$\dfrac{c(t)}{c_0}$")
    (l_phi,) = ax.plot(t, phi_ratio, label=r"$\dfrac{\varphi(t)}{\varphi_0}$")
    (l_psi,) = ax.plot(t, psi_ratio, label=r"$\dfrac{\psi(t)}{\psi_0}$")

    # Stage boundaries (only the meaningful ones)
    stage_lines = [0.0, rho_t1, load_t0, load_t1, t0, t1]
    for x in stage_lines:
        ax.axvline(
            x, color="black", linestyle="--", linewidth=1.6, alpha=0.75, zorder=0
        )

    # Stage labels
    spans = [
        (0.0, rho_t1, "Stage 1:\nGravity"),
        (load_t0, load_t1, "Stage 2:\nSurcharge"),
        (load_t1, t0, "Stage 3:\nHold"),
        (t0, t1, "Stage 4:\nSSR"),
    ]
    for a, b, txt in spans:
        ax.annotate(
            txt,
            xy=(0.5 * (a + b), 1.03),
            xycoords=("data", "axes fraction"),
            ha="center",
            va="bottom",
            fontsize=12,
            clip_on=False,
        )

    ax.set_xlabel(r"Pseudo-time $t$ / 1")
    ax.set_ylabel("Normalised value / 1")
    ax.set_xlim(0.0, float(t_end))
    ax.set_ylim(-0.03, 1.10)
    ax.grid(True, alpha=0.25)

    # Secondary axis: imposed trial factor
    ax2 = ax.twinx()
    (l_F,) = ax2.plot(
        t,
        F_plot,
        linestyle="--",
        linewidth=2.0,
        label=r"Imposed trial $F(t)=\dfrac{1}{s(t)}$",
    )
    ax2.set_ylabel(r"Trial reduction factor $F(t)$ / 1")
    ax2.set_ylim(0.95, max(2.3, 1.05 * float(np.nanmax(F_plot))))

    # Legends outside (no unused variables)
    handles_left = [l_rho, l_q, l_c, l_phi, l_psi]
    ax.legend(
        handles_left,
        [h.get_label() for h in handles_left],
        loc="center left",
        bbox_to_anchor=(1.18, 0.62),
        frameon=True,
    )
    ax2.legend(
        [l_F],
        [l_F.get_label()],
        loc="center left",
        bbox_to_anchor=(1.18, 0.20),
        frameon=True,
    )

    # Footer with initial values
    footer = (
        rf"$\rho_0={rho0:g}\,\mathrm{{kg/m^3}}$,  "
        rf"$q_0={q0:g}\,\mathrm{{Pa}}$,  "
        rf"$c_0={c0:g}\,\mathrm{{Pa}}$,  "
        rf"$\varphi_0={phi0_deg:g}^\circ$,  "
        rf"$\psi_0={psi0_deg:g}^\circ$,  "
        rf"$s_\min={smin:g}$"
    )
    fig.text(0.5, -0.01, footer, ha="center", va="bottom", fontsize=12)
    explanation = (
        "Curves show the staged loading and strength-reduction schedules versus pseudo-time. "
        "Vertical dashed lines mark stage transitions; the right axis shows the imposed trial "
        "factor F(t) used for SSR."
    )
    fig.text(0.5, -0.06, explanation, ha="center", va="bottom", fontsize=11)

    plt.subplots_adjust(top=0.86, right=0.78, bottom=0.18)

    png = out_dir / "ssr_staging_schedule.png"
    pdf = out_dir / "ssr_staging_schedule.pdf"
    fig.savefig(png, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")

    if show:
        plt.show()
    plt.close(fig)

    return png, pdf


# ============================================================
# OGS run
# ============================================================


def run_ogs(
    *,
    prj_path: Path,
    mesh_dir: Path,
    results_dir: Path,
    logfile: Path,
) -> None:
    prj = ot.Project(input_file=prj_path, output_file=prj_path)
    prj.run_model(args=f"-m {mesh_dir} -o {results_dir}", logfile=logfile)


# ============================================================
# Runner
# ============================================================
class MultiCaseRunner:
    def __init__(
        self, *, out_root: Path, tail_lines: int = 60, show_curves: bool = True
    ) -> None:
        self.out_root = Path(out_root)
        self.out_root.mkdir(parents=True, exist_ok=True)
        self.tail_lines = int(tail_lines)
        self.show_curves = bool(show_curves)

    def _tail(self, p: Path) -> str | None:
        try:
            if not p.exists():
                return None
            lines = p.read_text(errors="replace").splitlines()
            return "\n".join(lines[-self.tail_lines :]) if lines else None
        except Exception:
            return None

    def _check_mesh_folder(self, mesh_dir: Path, template_prj: Path) -> None:
        bulk = read_template_bulk_mesh_filename(template_prj)
        if not (mesh_dir / bulk).exists():
            msg = f"Bulk mesh '{bulk}' missing in mesh_dir: {mesh_dir}"
            raise FileNotFoundError(msg)

        expected = _all_vtu_mesh_names_in_prj(template_prj)
        missing = [fn for fn in sorted(expected) if not (mesh_dir / fn).exists()]
        if missing:
            missing_lines = "\n".join(f"  - {fn}" for fn in missing)
            msg = (
                "Mesh folder is missing VTU files required by template PRJ:\n"
                f"{missing_lines}\n\nmesh_dir = {mesh_dir}"
            )
            raise FileNotFoundError(msg)

    def run_one(
        self,
        *,
        basename: str,
        template_prj: Path,
        master: dict,
        case: CaseSpec,
    ) -> OGSResult:
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
        mesh_dir = (
            Path(case.mesh_dir) if case.mesh_dir is not None else (case_dir / "mesh")
        )
        results_dir = case_dir / "results"
        case_dir.mkdir(parents=True, exist_ok=True)
        results_dir.mkdir(parents=True, exist_ok=True)
        logfile = results_dir / "ogs.log"

        try:
            print(f"\n[RUN] {case_name}")
            print(f"      mesh_dir = {mesh_dir}")

            if not mesh_dir.exists():
                msg = f"mesh_dir does not exist: {mesh_dir}"
                raise FileNotFoundError(msg)

            self._check_mesh_folder(mesh_dir, template_prj)

            prj_out = case_dir / "project.prj"
            apply_cfg_to_prj(
                template_prj=template_prj,
                prj_out=prj_out,
                master=master,
                case=case,
            )

            plot_staging_and_ssr_curves(
                master=master,
                case_dir=case_dir,
                _case_name=case_name,
                show=self.show_curves,
            )

            run_ogs(
                prj_path=prj_out,
                mesh_dir=mesh_dir,
                results_dir=results_dir,
                logfile=logfile,
            )

            return OGSResult(
                name=case_name,
                case_dir=case_dir,
                ok=True,
                stage="ogs",
                logfile=logfile,
            )

        except Exception as e:
            (results_dir / "ogs_error_trace.txt").write_text(traceback.format_exc())
            return OGSResult(
                name=case_name,
                case_dir=case_dir,
                ok=False,
                stage="failed",
                error=f"{type(e).__name__}: {e}",
                logfile=logfile,
                log_tail=self._tail(logfile),
            )

    def run_all(
        self,
        *,
        basename: str,
        template_prj: Path,
        master: dict,
        cases: Sequence[CaseSpec],
    ) -> list[OGSResult]:
        out: list[OGSResult] = []
        for c in cases:
            r = self.run_one(
                basename=basename, template_prj=template_prj, master=master, case=c
            )
            print("[OK]" if r.ok else "[FAILED]", r.name)

            if not r.ok:
                print("  error:", r.error)
                trace = r.case_dir / "results" / "ogs_error_trace.txt"
                if trace.exists():
                    print("  trace:", trace)

                if r.log_tail:
                    print("  --- ogs.log tail ---")
                    print(r.log_tail)
                    print("  --------------------")

            out.append(r)
        return out
