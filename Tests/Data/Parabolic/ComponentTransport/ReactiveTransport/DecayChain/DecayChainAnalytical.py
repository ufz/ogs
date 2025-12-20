# SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
# SPDX-License-Identifier: BSD-3-Clause

import types

import numpy as np
from scipy import special


# cf. https://stackoverflow.com/a/31047259 and https://stackoverflow.com/a/54249582
def noglobals(f):
    return types.FunctionType(
        f.__code__, globals().copy(), f.__name__, f.__defaults__, f.__closure__
    )


def computeProduct(j, i, k, c_inlet):
    value = 1
    for l in range(j, i):  # noqa: E741
        value *= k[l] / (k[l] - k[i]) * c_inlet[j]

    return value


def computeInitialAuxiliaryVariable(c_inlet, k):
    a_inlet = np.empty(0)

    for i in range(len(c_inlet)):
        value = c_inlet[i]
        if i > 0:
            for j in range(i):
                value += computeProduct(j, i, k, c_inlet)
        a_inlet = np.append(a_inlet, value)

    return a_inlet


def computeAnalyticalSolution(x, t, c_0, k, v, D):
    t *= 3.1536e7  # unit conversion from year to second

    beta = (v**2 / 4 / D**2 + k / D) ** 0.5
    return (
        c_0
        / 2
        * np.exp(v * x / 2 / D)
        * (
            np.exp(-beta * x)
            * special.erfc((x - (v**2 + 4 * k * D) ** 0.5 * t) / 2 / (D * t) ** 0.5)
            + np.exp(beta * x)
            * special.erfc((x + (v**2 + 4 * k * D) ** 0.5 * t) / 2 / (D * t) ** 0.5)
        )
    )


def computeGradAnalyticalSolution(x, t, c_0, k, v, D):
    t *= 3.1536e7  # unit conversion from year to second

    assert v == 0  # v != 0 contributions not yet implemented

    beta = (v**2 / 4 / D**2 + k / D) ** 0.5

    sdt = np.sqrt(D * t)
    embx = np.exp(-beta * x)
    epbx = np.exp(+beta * x)
    argm = (x - (v**2 + 4 * k * D) ** 0.5 * t) / (2 * sdt)
    argp = (x + (v**2 + 4 * k * D) ** 0.5 * t) / (2 * sdt)
    erfcm = special.erfc(argm)
    erfcp = special.erfc(argp)
    expm = np.exp(-(argm**2))
    expp = np.exp(-(argp**2))

    return (
        0.5
        * c_0
        * (
            -beta * embx * erfcm
            - embx * 2 / np.sqrt(np.pi) * expm / (2 * sdt)
            + beta * epbx * erfcp
            - epbx * 2 / np.sqrt(np.pi) * expp / (2 * sdt)
        )
    )


@noglobals
def computeConcentrations(x, t, v, D, k, radionuclides, c_inlet):
    c = {}
    a = {}

    a_inlet = dict(
        zip(
            radionuclides,
            computeInitialAuxiliaryVariable(c_inlet, list(k.values())),
            strict=True,
        )
    )

    c["[Cm-247]"] = computeAnalyticalSolution(
        x, t, a_inlet["[Cm-247]"], k["[Cm-247]"], v, D
    )

    a["[Am-243]"] = computeAnalyticalSolution(
        x, t, a_inlet["[Am-243]"], k["[Am-243]"], v, D
    )
    c["[Am-243]"] = (
        a["[Am-243]"] - k["[Cm-247]"] / (k["[Cm-247]"] - k["[Am-243]"]) * c["[Cm-247]"]
    )

    a["[Pu-239]"] = computeAnalyticalSolution(
        x, t, a_inlet["[Pu-239]"], k["[Pu-239]"], v, D
    )
    c["[Pu-239]"] = (
        a["[Pu-239]"]
        - k["[Cm-247]"]
        / (k["[Cm-247]"] - k["[Pu-239]"])
        * k["[Am-243]"]
        / (k["[Am-243]"] - k["[Pu-239]"])
        * c["[Cm-247]"]
        - k["[Am-243]"] / (k["[Am-243]"] - k["[Pu-239]"]) * c["[Am-243]"]
    )

    a["[U-235]"] = computeAnalyticalSolution(
        x, t, a_inlet["[U-235]"], k["[U-235]"], v, D
    )
    c["[U-235]"] = (
        a["[U-235]"]
        - k["[Cm-247]"]
        / (k["[Cm-247]"] - k["[U-235]"])
        * k["[Am-243]"]
        / (k["[Am-243]"] - k["[U-235]"])
        * k["[Pu-239]"]
        / (k["[Pu-239]"] - k["[U-235]"])
        * c["[Cm-247]"]
        - k["[Am-243]"]
        / (k["[Am-243]"] - k["[U-235]"])
        * k["[Pu-239]"]
        / (k["[Pu-239]"] - k["[U-235]"])
        * c["[Am-243]"]
        - k["[Pu-239]"] / (k["[Pu-239]"] - k["[U-235]"]) * c["[Pu-239]"]
    )

    a["[Pa-231]"] = computeAnalyticalSolution(
        x, t, a_inlet["[Pa-231]"], k["[Pa-231]"], v, D
    )
    c["[Pa-231]"] = (
        a["[Pa-231]"]
        - k["[Cm-247]"]
        / (k["[Cm-247]"] - k["[Pa-231]"])
        * k["[Am-243]"]
        / (k["[Am-243]"] - k["[Pa-231]"])
        * k["[Pu-239]"]
        / (k["[Pu-239]"] - k["[Pa-231]"])
        * k["[U-235]"]
        / (k["[U-235]"] - k["[Pa-231]"])
        * c["[Cm-247]"]
        - k["[Am-243]"]
        / (k["[Am-243]"] - k["[Pa-231]"])
        * k["[Pu-239]"]
        / (k["[Pu-239]"] - k["[Pa-231]"])
        * k["[U-235]"]
        / (k["[U-235]"] - k["[Pa-231]"])
        * c["[Am-243]"]
        - k["[Pu-239]"]
        / (k["[Pu-239]"] - k["[Pa-231]"])
        * k["[U-235]"]
        / (k["[U-235]"] - k["[Pa-231]"])
        * c["[Pu-239]"]
        - k["[U-235]"] / (k["[U-235]"] - k["[Pa-231]"]) * c["[U-235]"]
    )

    a["[Ac-227]"] = computeAnalyticalSolution(
        x, t, a_inlet["[Ac-227]"], k["[Ac-227]"], v, D
    )
    c["[Ac-227]"] = (
        a["[Ac-227]"]
        - k["[Cm-247]"]
        / (k["[Cm-247]"] - k["[Ac-227]"])
        * k["[Am-243]"]
        / (k["[Am-243]"] - k["[Ac-227]"])
        * k["[Pu-239]"]
        / (k["[Pu-239]"] - k["[Ac-227]"])
        * k["[U-235]"]
        / (k["[U-235]"] - k["[Ac-227]"])
        * k["[Pa-231]"]
        / (k["[Pa-231]"] - k["[Ac-227]"])
        * c["[Cm-247]"]
        - k["[Am-243]"]
        / (k["[Am-243]"] - k["[Ac-227]"])
        * k["[Pu-239]"]
        / (k["[Pu-239]"] - k["[Ac-227]"])
        * k["[U-235]"]
        / (k["[U-235]"] - k["[Ac-227]"])
        * k["[Pa-231]"]
        / (k["[Pa-231]"] - k["[Ac-227]"])
        * c["[Am-243]"]
        - k["[Pu-239]"]
        / (k["[Pu-239]"] - k["[Ac-227]"])
        * k["[U-235]"]
        / (k["[U-235]"] - k["[Ac-227]"])
        * k["[Pa-231]"]
        / (k["[Pa-231]"] - k["[Ac-227]"])
        * c["[Pu-239]"]
        - k["[U-235]"]
        / (k["[U-235]"] - k["[Ac-227]"])
        * k["[Pa-231]"]
        / (k["[Pa-231]"] - k["[Ac-227]"])
        * c["[U-235]"]
        - k["[Pa-231]"] / (k["[Pa-231]"] - k["[Ac-227]"]) * c["[Pa-231]"]
    )

    return c


@noglobals
def computeGradients(x, t, v, D, k, radionuclides, c_inlet):
    ###Analytical solution###
    grad_c = {}
    grad_a = {}

    a_inlet = dict(
        zip(
            radionuclides,
            computeInitialAuxiliaryVariable(c_inlet, list(k.values())),
            strict=True,
        )
    )

    grad_c["[Cm-247]"] = computeGradAnalyticalSolution(
        x, t, a_inlet["[Cm-247]"], k["[Cm-247]"], v, D
    )

    grad_a["[Am-243]"] = computeGradAnalyticalSolution(
        x, t, a_inlet["[Am-243]"], k["[Am-243]"], v, D
    )
    grad_c["[Am-243]"] = (
        grad_a["[Am-243]"]
        - k["[Cm-247]"] / (k["[Cm-247]"] - k["[Am-243]"]) * grad_c["[Cm-247]"]
    )

    grad_a["[Pu-239]"] = computeGradAnalyticalSolution(
        x, t, a_inlet["[Pu-239]"], k["[Pu-239]"], v, D
    )
    grad_c["[Pu-239]"] = (
        grad_a["[Pu-239]"]
        - k["[Cm-247]"]
        / (k["[Cm-247]"] - k["[Pu-239]"])
        * k["[Am-243]"]
        / (k["[Am-243]"] - k["[Pu-239]"])
        * grad_c["[Cm-247]"]
        - k["[Am-243]"] / (k["[Am-243]"] - k["[Pu-239]"]) * grad_c["[Am-243]"]
    )

    grad_a["[U-235]"] = computeGradAnalyticalSolution(
        x, t, a_inlet["[U-235]"], k["[U-235]"], v, D
    )
    grad_c["[U-235]"] = (
        grad_a["[U-235]"]
        - k["[Cm-247]"]
        / (k["[Cm-247]"] - k["[U-235]"])
        * k["[Am-243]"]
        / (k["[Am-243]"] - k["[U-235]"])
        * k["[Pu-239]"]
        / (k["[Pu-239]"] - k["[U-235]"])
        * grad_c["[Cm-247]"]
        - k["[Am-243]"]
        / (k["[Am-243]"] - k["[U-235]"])
        * k["[Pu-239]"]
        / (k["[Pu-239]"] - k["[U-235]"])
        * grad_c["[Am-243]"]
        - k["[Pu-239]"] / (k["[Pu-239]"] - k["[U-235]"]) * grad_c["[Pu-239]"]
    )

    grad_a["[Pa-231]"] = computeGradAnalyticalSolution(
        x, t, a_inlet["[Pa-231]"], k["[Pa-231]"], v, D
    )
    grad_c["[Pa-231]"] = (
        grad_a["[Pa-231]"]
        - k["[Cm-247]"]
        / (k["[Cm-247]"] - k["[Pa-231]"])
        * k["[Am-243]"]
        / (k["[Am-243]"] - k["[Pa-231]"])
        * k["[Pu-239]"]
        / (k["[Pu-239]"] - k["[Pa-231]"])
        * k["[U-235]"]
        / (k["[U-235]"] - k["[Pa-231]"])
        * grad_c["[Cm-247]"]
        - k["[Am-243]"]
        / (k["[Am-243]"] - k["[Pa-231]"])
        * k["[Pu-239]"]
        / (k["[Pu-239]"] - k["[Pa-231]"])
        * k["[U-235]"]
        / (k["[U-235]"] - k["[Pa-231]"])
        * grad_c["[Am-243]"]
        - k["[Pu-239]"]
        / (k["[Pu-239]"] - k["[Pa-231]"])
        * k["[U-235]"]
        / (k["[U-235]"] - k["[Pa-231]"])
        * grad_c["[Pu-239]"]
        - k["[U-235]"] / (k["[U-235]"] - k["[Pa-231]"]) * grad_c["[U-235]"]
    )

    grad_a["[Ac-227]"] = computeGradAnalyticalSolution(
        x, t, a_inlet["[Ac-227]"], k["[Ac-227]"], v, D
    )
    grad_c["[Ac-227]"] = (
        grad_a["[Ac-227]"]
        - k["[Cm-247]"]
        / (k["[Cm-247]"] - k["[Ac-227]"])
        * k["[Am-243]"]
        / (k["[Am-243]"] - k["[Ac-227]"])
        * k["[Pu-239]"]
        / (k["[Pu-239]"] - k["[Ac-227]"])
        * k["[U-235]"]
        / (k["[U-235]"] - k["[Ac-227]"])
        * k["[Pa-231]"]
        / (k["[Pa-231]"] - k["[Ac-227]"])
        * grad_c["[Cm-247]"]
        - k["[Am-243]"]
        / (k["[Am-243]"] - k["[Ac-227]"])
        * k["[Pu-239]"]
        / (k["[Pu-239]"] - k["[Ac-227]"])
        * k["[U-235]"]
        / (k["[U-235]"] - k["[Ac-227]"])
        * k["[Pa-231]"]
        / (k["[Pa-231]"] - k["[Ac-227]"])
        * grad_c["[Am-243]"]
        - k["[Pu-239]"]
        / (k["[Pu-239]"] - k["[Ac-227]"])
        * k["[U-235]"]
        / (k["[U-235]"] - k["[Ac-227]"])
        * k["[Pa-231]"]
        / (k["[Pa-231]"] - k["[Ac-227]"])
        * grad_c["[Pu-239]"]
        - k["[U-235]"]
        / (k["[U-235]"] - k["[Ac-227]"])
        * k["[Pa-231]"]
        / (k["[Pa-231]"] - k["[Ac-227]"])
        * grad_c["[U-235]"]
        - k["[Pa-231]"] / (k["[Pa-231]"] - k["[Ac-227]"]) * grad_c["[Pa-231]"]
    )

    return grad_c
