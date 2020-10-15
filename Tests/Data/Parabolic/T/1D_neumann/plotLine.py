#!/usr/bin/env python

import temperature_analytical
from matplotlib.legend import Legend
from vtk import *
from vtk.numpy_interface import dataset_adapter as dsa
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import scipy.optimize as sp


def probeFileAlongLine(filename):
    reader = vtkXMLUnstructuredGridReader()
    reader.SetFileName(filename)

    line = vtkLineSource()
    line.SetResolution(10000)
    line.SetPoint1(0, 0, 0)
    line.SetPoint2(60, 0, 0)
    line_probe = vtkProbeFilter()
    line_probe.SetInputConnection(line.GetOutputPort())
    line_probe.SetSourceConnection(reader.GetOutputPort())
    line_probe.Update()

    return dsa.WrapDataObject(line_probe.GetOutput())


def convertPointDataToPandasDF(line_data):
    point_data = line_data.GetPointData()
    data = {}
    for k, v in zip(point_data.keys(), point_data.values()):
        data[k] = vtk_to_numpy(v)

    coords = line_data.GetPoints()[:, 0]  # x-coordinates
    df = pd.concat(
        [pd.DataFrame(v) for k, v in data.items()], axis=1, keys=list(data.keys())
    )
    df.index = coords
    return df


def setupMatplotlib(height=8.0, width=6.0):
    plt.rcParams["xtick.direction"] = "out"
    plt.rcParams["ytick.direction"] = "out"
    plt.rcParams["lines.linewidth"] = 2.0
    plt.rcParams["lines.color"] = "black"
    plt.rcParams["legend.frameon"] = True
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["legend.fontsize"] = 8
    plt.rcParams["font.size"] = 12
    # For ipython notebook display set default values.
    # plt.rcParams['lines.markersize'] = 12
    plt.rcParams["figure.figsize"] = (height, width)
    plt.rcParams["grid.linewidth"] = 1

    # General settings used by display and print contexts.
    plt.rcParams["axes.axisbelow"] = True
    grid_line_color = "0.95"
    plt.rcParams["grid.color"] = grid_line_color
    plt.rcParams["grid.linestyle"] = "-"
    plt.rc("text", usetex=True)
    plt.rc("font", family="serif")


def commonFormat(ax_el, centerx=None, centery=None):
    ax_el.grid(True)
    ax_el.spines["top"].set_visible(0)
    ax_el.spines["right"].set_visible(0)
    ax_el.spines["bottom"].set_linewidth(0.5)
    ax_el.spines["left"].set_linewidth(0.5)
    ax_el.xaxis.set_ticks_position("bottom")
    ax_el.yaxis.set_ticks_position("left")
    if (centerx is not None) and (centery is not None):
        ax_el.spines["left"].set_position(("data", centerx))
        ax_el.spines["bottom"].set_position(("data", centery))
        ax_el.spines["right"].set_position(("data", centerx - 1))
        ax_el.spines["top"].set_position(("data", centery - 1))


def plotWithFormatting(ax, data, formatting, with_labels=True):
    values = data[formatting["column"]]
    if "component" in formatting:
        values = values[formatting["component"]]

    # print("plot", formatting['column'], values)
    return ax.plot(
        data.index,
        values,
        ls=formatting["ls"],
        color=formatting["color"],
        lw=1.0,
        label=formatting["label"] if with_labels else "",
    )


def plotTemperature(ax, data, with_labels=True):
    formattings = [
        {
            "column": ("newton", "temperature"),
            "component": 0,
            "color": "green",
            "ls": "-",
            "label": "no mass lumping",
        },
        {
            "column": ("newton_masslumping", "temperature"),
            "component": 0,
            "color": "blue",
            "ls": "-",
            "label": "with mass lumping",
        },
        {
            "column": ("newton", "analytical"),
            "color": "#e41a1c",
            "ls": "-",
            "label": "analytical",
        },
    ]
    lines = []
    for formatting in formattings:
        lines += plotWithFormatting(ax, data, formatting, with_labels)
    return lines


def plotErrors(ax, data, with_labels=True):
    formattings = [
        {
            "column": ("newton", "error"),
            "color": "green",
            "ls": "-",
            "label": "$e_\mathrm{newton}$",
        },
        {
            "column": ("newton_masslumping", "error"),
            "color": "blue",
            "ls": "-",
            "label": "$e_\mathrm{newton}^\mathrm{ML}$",
        },
    ]
    lines = []
    for formatting in formattings:
        lines += plotWithFormatting(ax, data, formatting, with_labels)
    return lines


def plotCases(data, output):
    d = data[data[("picard_masslumping", "error")].abs() > 1e-6]
    setupMatplotlib(6, 7)
    plt.rcParams["legend.fontsize"] = 10
    plt.rcParams["font.size"] = 14
    fig, axes = plt.subplots(nrows=2, ncols=1)
    lines = plotTemperature(axes[0], d, with_labels=True)
    axes[0].set_ylabel("Temperature / K", fontsize=18)
    axes[0].legend(loc="best", fontsize=12, ncol=1, bbox_to_anchor=(0.5, 0.5, 0.5, 0.5))

    lines += plotErrors(axes[1], d, with_labels=False)
    axes[1].set_ylabel("Error / K", fontsize=18)
    axes[1].set_xlabel("$x$ / m", fontsize=16)

    for ax in axes:
        commonFormat(ax)
        # ax.set_xlim(left=0.)
        # ax.set_xlim(right=5.)
        fig.tight_layout()
    fig.savefig(output + ".png", dpi=150)
    plt.close("all")
    return None


def plotPicardTemperature(ax, data, with_labels=True):
    formattings = [
        {
            "column": ("picard", "temperature"),
            "component": 0,
            "color": "blue",
            "ls": "-",
            "label": "picard, no mass lumping",
        },
        {
            "column": ("picard_masslumping", "temperature"),
            "component": 0,
            "color": "green",
            "ls": "-",
            "label": "picard, with mass lumping",
        },
    ]
    lines = []
    for formatting in formattings:
        lines += plotWithFormatting(ax, data, formatting, with_labels)
    return lines


def plotNewtonTemperature(ax, data, with_labels=True):
    formattings = [
        {
            "column": ("newton", "temperature"),
            "component": 0,
            "color": "orange",
            "ls": "-.",
            "label": "newton, no mass lumping",
        },
        {
            "column": ("newton_masslumping", "temperature"),
            "component": 0,
            "color": "red",
            "ls": "-.",
            "label": "newton, with mass lumping",
        },
    ]
    lines = []
    for formatting in formattings:
        lines += plotWithFormatting(ax, data, formatting, with_labels)
    return lines


def plotNewtonVsPicardErrors(ax, data, with_labels=True):
    formattings = [
        {"column": "no_ml", "color": "green", "ls": "-", "label": "no mass lumping"},
        {"column": "ml", "color": "blue", "ls": "-", "label": "with mass lumping"},
    ]
    lines = []
    for formatting in formattings:
        lines += plotWithFormatting(ax, data, formatting, with_labels)
    return lines


def plotNewtonVsPicard(data, output):
    data["no_ml"] = (
        data[("newton", "temperature", 0)] - data[("picard", "temperature", 0)]
    )
    data["ml"] = (
        data[("newton_masslumping", "temperature", 0)]
        - data[("picard_masslumping", "temperature", 0)]
    )
    # data['select'] = (data['ml'].abs() > 1e-16) | (data['no_ml'].abs() > 1e-16)
    # print(d[d].index[-1])

    setupMatplotlib(6, 7)
    plt.rcParams["legend.fontsize"] = 10
    plt.rcParams["font.size"] = 14
    fig, axes = plt.subplots(nrows=2, ncols=1)
    lines = plotPicardTemperature(axes[0], data, with_labels=True)
    lines = plotNewtonTemperature(axes[0], data, with_labels=True)
    axes[0].set_ylabel("Temperature / K", fontsize=18)
    axes[0].legend(loc="best", fontsize=12, ncol=1, bbox_to_anchor=(0.5, 0.5, 0.5, 0.5))

    lines += plotNewtonVsPicardErrors(axes[1], data, with_labels=True)
    axes[1].legend(loc="best", fontsize=12, ncol=1, bbox_to_anchor=(0.5, 0.5, 0.5, 0.5))
    axes[1].set_ylabel("$T_\mathrm{diff}$ (newton - picard) / K", fontsize=18)
    axes[1].set_xlabel("$x$ / m", fontsize=16)

    for ax in axes:
        commonFormat(ax)
        # ax.set_xlim(left=0.)
        # ax.set_xlim(right=5.)
        fig.tight_layout()
    fig.savefig(output + ".png", dpi=300)
    plt.close("all")
    return None


def singleTimeStep(ts):
    dt = 78125
    t = ts * dt
    cases = ["picard", "newton", "picard_masslumping", "newton_masslumping"]
    data = []
    for c in cases:
        line = probeFileAlongLine(c + "_ts_" + str(ts) + "_t_" + str(t) + ".000000.vtu")
        df = convertPointDataToPandasDF(line)
        df["analytical"] = temperature_analytical.temperature(df.index, t)
        df["error"] = df["temperature"].sub(df["analytical"], axis=0)
        print(
            "Error for ts",
            ts,
            "case",
            c,
            "is",
            np.sqrt((df["error"] ** 2).sum() * 60 / 10000),
        )
        data.append(df[["temperature", "error", "analytical"]])

    df = pd.concat(data, keys=cases, names=["case"], axis=1)
    plotCases(
        df[df[("picard_masslumping", "error")].abs() > 1e-6],
        "temperature_error_ts_" + str(ts) + "_t_" + str(t),
    )
    plotNewtonVsPicard(df, "picard_vs_newton_ts_" + str(ts) + "_t_" + str(t))


def main():
    for ts in [1, 3, 65, 405]:
        singleTimeStep(ts)
        print()


if __name__ == "__main__":
    main()
