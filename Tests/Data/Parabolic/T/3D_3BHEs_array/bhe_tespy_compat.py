# SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
# SPDX-License-Identifier: BSD-3-Clause

"""Shared TESPy glue for the BHE pipe-network benchmarks.

``install_workarounds()`` has to be called *before* the first ``import tespy``,
because TESPy pulls in CoolProp at module import time.

An identical copy of this file lives next to every benchmark that needs it, so
that a benchmark directory stays self-contained when a user copies it as a
template.  Keep the copies byte-identical (``cmp`` is a valid drift check).

Tested against TESPy 0.9.16 .. 0.11.x.
"""

import sys
import types

# Module-level names that TESPy's CoolProp back end looks up while importing.
# The values are never used, because all networks here are configured for the
# IAPWSWrapper fluid property engine.
_COOLPROP_CONSTANTS = (
    "iT_min",
    "iT_max",
    "iP_min",
    "iP_max",
    "iP_critical",
    "iT_critical",
    "iT_freeze",
    "imolar_mass",
    "HmassP_INPUTS",
    "PSmass_INPUTS",
    "PT_INPUTS",
    "QT_INPUTS",
    "PQ_INPUTS",
    "iphase_twophase",
    "iphase_liquid",
    "iphase_gas",
    "iphase_supercritical_gas",
)


def install_workarounds():
    """Make TESPy usable inside OGS's embedded Python interpreter.

    Call this before importing anything from ``tespy``.  Calling it more than
    once is harmless.
    """
    # Order matters: the CoolProp mock has to be in place before the patch
    # below, because it imports TESPy.
    _mock_coolprop()
    _patch_iapws_viscosity()


def _mock_coolprop():
    """Install a stub CoolProp module to prevent a segfault inside OGS.

    CoolProp 7.2.0 crashes in ``PredefinedMixturesLibrary::load_from_JSON()``
    when loaded inside OGS's embedded interpreter.  TESPy imports CoolProp
    unconditionally at module level, so the only way to avoid the crash is to
    make the import find a stub instead.  The networks use IAPWSWrapper (see
    the ``engine`` entries in ``pre/tespy_nw*.json``), so nothing in the stub
    is ever called.
    """
    mock = types.ModuleType("CoolProp")
    mock.AbstractState = type("AbstractState", (), {})
    # TESPy imports HAPropsSI at module level in
    # tespy/tools/fluid_properties/functions.py; a dummy satisfies the import.
    mock.HAPropsSI = lambda *args, **kwargs: 0.0
    for name in _COOLPROP_CONSTANTS:
        setattr(mock, name, 0)
    sys.modules["CoolProp"] = mock
    sys.modules["CoolProp.CoolProp"] = mock


def _patch_iapws_viscosity():
    """Give IAPWSWrapper a viscosity fallback in the two-phase region.

    ``iapws.IAPWS97(...).mu`` returns ``None`` for two-phase states.  The
    benchmark loops carry liquid water at 283..305 K, so two-phase states are
    not physical here -- they only show up as transient excursions while the
    TESPy solver is still far from the solution.  Falling back to the
    saturated liquid viscosity keeps the iteration going and does not affect
    the converged result.

    The fallback can itself return ``None`` if the saturation lookup fails; in
    that case TESPy raises further down, which is the intended behaviour.
    """
    # imported here, not at module level, so that _mock_coolprop() gets to run
    # before anything pulls in TESPy
    from tespy.tools.fluid_properties import wrappers  # noqa: PLC0415

    IAPWSWrapper = wrappers.IAPWSWrapper

    if not hasattr(IAPWSWrapper, "viscosity_ph"):
        msg = (
            "tespy.tools.fluid_properties.wrappers.IAPWSWrapper.viscosity_ph "
            "no longer exists, so the two-phase viscosity workaround cannot "
            "be applied. Check whether iapws now returns a viscosity in the "
            "two-phase region and drop this patch."
        )
        raise RuntimeError(msg)

    def viscosity_ph(self, p, h):
        state = self.AS(P=p / 1e6, h=h / 1e3)
        if state.mu is not None:
            return state.mu
        return self.AS(P=p / 1e6, x=0).mu

    IAPWSWrapper.viscosity_ph = viscosity_ph


def find_network_inlet(nw):
    """Return the connection that feeds the network.

    That is the connection leaving the ``Source`` (open loop) or the
    ``CycleCloser`` (closed loop).  Looking the component up by type rather
    than by a hard-coded label keeps this working when the network is renamed
    in ``pre/3bhes*.py``, and fails loudly instead of silently binding
    nothing.
    """
    # imported here, not at module level, so that _mock_coolprop() gets to run
    # before anything pulls in TESPy
    from tespy.components import CycleCloser, Source  # noqa: PLC0415

    matches = [
        c for c in nw.conns["object"] if isinstance(c.source, (Source, CycleCloser))
    ]
    if len(matches) != 1:
        labels = sorted(c.label for c in matches)
        msg = (
            f"Expected exactly one network inlet connection (leaving a Source "
            f"or a CycleCloser), found {len(matches)}: {labels}."
        )
        raise RuntimeError(msg)
    return matches[0]


def get_component(nw, label):
    """Return the component with ``label``, raising if it does not exist.

    ``Network.get_comp()`` returns ``None`` for an unknown label (a
    ``FutureWarning`` since TESPy 0.11), which turns a typo in the user
    settings into an ``AttributeError`` much later on.
    """
    comp = nw.get_comp(label)
    if comp is None:
        known = sorted(nw.comps.index)
        msg = f"No component labelled {label!r} in the network. Known: {known}."
        raise RuntimeError(msg)
    return comp
