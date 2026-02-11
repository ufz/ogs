from __future__ import annotations

import math


def _clamp(x: float, lo: float, hi: float) -> float:
    return max(lo, min(hi, x))


def _require(cond: bool, msg: str) -> None:
    if not cond:
        err = msg
        raise ValueError(err)


# =========================
# Schedule parameters
# =========================
def stage_params(MASTER: dict) -> tuple[float, float, float]:
    stages = MASTER.get("stages")
    if stages is None:
        msg = "MASTER must contain 'stages' for rho/load expressions."
        raise KeyError(msg)

    rho_t1 = float(stages["rho_t1"])
    load_t0 = float(stages["load_t0"])
    load_t1 = float(stages["load_t1"])

    msg = f"Invalid stages: rho_t1 must be > 0, got {rho_t1}."
    _require(rho_t1 > 0.0, msg)

    msg = f"Invalid stages: load_t1 must be > load_t0, got {load_t0}..{load_t1}."
    _require(load_t1 > load_t0, msg)

    return rho_t1, load_t0, load_t1


def ssr_params(MASTER: dict) -> tuple[float, float, float, float]:
    ssr = MASTER.get("ssr")
    if ssr is None:
        msg = "MASTER must contain 'ssr' with keys: t0, t1, smin."
        raise KeyError(msg)

    t0 = float(ssr["t0"])
    t1 = float(ssr["t1"])
    smin = float(ssr["smin"])

    msg = f"Invalid SSR: t1 must be > t0, got {t0}..{t1}."
    _require(t1 > t0, msg)

    msg = f"Invalid SSR: smin must be in (0,1], got {smin}."
    _require(0.0 < smin <= 1.0, msg)

    return t0, t1, smin, t1 - t0


# =========================
# Naming helpers
# =========================
def _float_key(x: float) -> str:
    return f"{x:.12g}".replace(".", "p").replace("-", "m")


def _sci_key(x: float) -> str:
    s = f"{x:.0e}" if x != 0 else "0"
    return s.replace("+", "")


def case_folder_name(
    basename: str,
    *,
    element_order: int,
    h: float,
    dtmin: float,
    tol: float,
    load_tag: str,
    label: str,
) -> str:
    ord_tag = "lin" if element_order == 1 else "quad"
    return (
        f"{basename}"
        f"__ord-{ord_tag}"
        f"__h-{_float_key(h)}"
        f"__dtmin-{_sci_key(dtmin)}"
        f"__tol-{_sci_key(tol)}"
        f"__load-{load_tag}"
        f"__{label}"
    )


# =========================
# SSR schedule
# =========================
def _ssr_strength_multiplier_linear(
    t: float, *, t0: float, t1: float, smin: float
) -> float:
    """
    Strength multiplier s(t)=1/F(t):
      - s=1 for t<=t0
      - linearly decreases to smin at t=t1
      - clamped to [smin, 1]
    """
    msg = f"Invalid SSR: t1 must be > t0, got {t0}..{t1}."
    _require(t1 > t0, msg)

    msg = f"Invalid SSR: smin must be in (0,1], got {smin}."
    _require(0.0 < smin <= 1.0, msg)

    if t <= t0:
        return 1.0
    if t >= t1:
        return float(smin)
    s = 1.0 - (1.0 - float(smin)) * (t - t0) / (t1 - t0)
    return _clamp(s, float(smin), 1.0)


def ssr_schedule(t: float, MASTER: dict) -> tuple[float, float, float, float, float]:
    """
    Returns:
      F(t), s(t), c(t), phi_deg(t), psi_deg(t)

    Uses:
      MASTER["ssr"] = {"t0":..., "t1":..., "smin":...}
      MASTER["material"]["parameters"]["scalar"] = {...}
    """
    ssr_t0, ssr_t1, ssr_smin, _dt = ssr_params(MASTER)

    t = float(t)
    s = _ssr_strength_multiplier_linear(t, t0=ssr_t0, t1=ssr_t1, smin=ssr_smin)
    F = 1.0 / s

    scalars = MASTER["material"]["parameters"]["scalar"]
    c0 = float(scalars["cohesion_rs0"])
    phi0_deg = float(scalars["FrictionAngle_rs0"])
    psi0_deg = float(scalars.get("DilatancyAngle_rs0", 0.0))

    phi0 = math.radians(phi0_deg)
    psi0 = math.radians(psi0_deg)

    c = c0 * s
    phi_deg = math.degrees(math.atan(math.tan(phi0) * s))
    psi_deg = math.degrees(math.atan(math.tan(psi0) * s))

    return F, s, c, phi_deg, psi_deg


# =================================
# Expression generators for the PRJ
# =================================
def ssr_strength_expr(MASTER: dict) -> str:
    """
    Returns expression string for s(t) (strength multiplier), clamped.
    """
    t0, t1, smin, _dt = ssr_params(MASTER)

    return (
        f"max({smin:.16g}, min(1.0, "
        f"1.0 - (1.0-{smin:.16g}) * (t-({t0:.16g})) / (({t1:.16g})-({t0:.16g}))"
        f"))"
    )


def rho_scale_expr(MASTER: dict) -> str:
    """
    rho_scale(t) = clamp(t/rho_t1, 0, 1)
    """
    t1, _load_t0, _load_t1 = stage_params(MASTER)
    return f"max(0.0, min(1.0, t/({t1:.16g})))"


def load_scale_expr(MASTER: dict) -> str:
    """
    load_scale(t) = clamp((t-load_t0)/(load_t1-load_t0), 0, 1)
    """
    _rho_t1, t0, t1 = stage_params(MASTER)
    return f"max(0.0, min(1.0, (t-({t0:.16g})) / (({t1:.16g})-({t0:.16g}))))"


def rho_sr_expr(MASTER: dict, *, rho0: float) -> str:
    """rho_sr = rho0 * rho_scale(t)"""
    return f"{rho0:.16g} * ({rho_scale_expr(MASTER)})"


def top_load_expr(MASTER: dict, *, top0: float) -> str:
    """top_load = top0 * load_scale(t)"""
    return f"{top0:.16g} * ({load_scale_expr(MASTER)})"
