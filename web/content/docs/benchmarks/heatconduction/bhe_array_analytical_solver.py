###
# Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
# Distributed under a Modified BSD License.
# See accompanying file LICENSE.txt or
# http://www.opengeosys.org/project/license
###
"""
Python script for superposition principle based
analytical solution to the temporal temperature
change at an arbitrary location during the operation
of a BHE array system proposed by Bayer2014

Bayer, P., de Paly, M., & Beck, M. (2014).
Strategic optimization of borehole heat exchanger field
for seasonal geothermal heating and cooling.
Applied Energy, 136, 445–453.
https://doi.org/10.1016/j.apenergy.2014.09.029

Author: Shuang Chen
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy import special as sp
import math

#%% input parameters
# source term coordinates
po_x = np.array(
    [
        40,
        40,
        40,
        40,
        40,
        45,
        45,
        45,
        45,
        45,
        50,
        50,
        50,
        50,
        50,
        55,
        55,
        55,
        55,
        55,
        60,
        60,
        60,
        60,
        60,
    ],
    dtype=float,
).reshape(-1, 1)
po_y = np.array(
    [
        60,
        55,
        50,
        45,
        40,
        60,
        55,
        50,
        45,
        40,
        60,
        55,
        50,
        45,
        40,
        60,
        55,
        50,
        45,
        40,
        60,
        55,
        50,
        45,
        40,
    ],
    dtype=float,
).reshape(-1, 1)

po_dist_to_referencepo = np.zeros([25, 1])
Temp_po_to_referencepo = np.zeros([25, 1])

# specific heat exchange rate in Watt/metre in each month:
q1 = -35
q2 = -35
q3 = -35
q4 = -35
q5 = 0
q6 = 0
q7 = 0
q8 = 0
q9 = 0
q10 = 0
q11 = 0
q12 = 0

time_trans = 30 * 24 * 60 * 60  # time stepping val

# parameters
T0 = 10  # initial subsurface temperature in degreeC

lamda_sp = 2.0  # solid thermal conductivity
density_sp = 1500  # solid density
cp_sp = 1950  # solid heat capacity

alpha_sp = lamda_sp / (density_sp * cp_sp)  # thermal dispersion tensor

# reference plot line A-A'section coordinates.
# A-A'section located 1 m away from the diagonal.
point_x = np.arange(0.5, 100.1, 0.5)
point_y = np.arange(0, 99.6, 0.5)

numbhe = 25  # BHE number in the array
num_refer_points = len(point_x)

# thermal load curve
qq = np.array(
    [
        q1,
        q2,
        q3,
        q4,
        q5,
        q6,
        q7,
        q8,
        q9,
        q10,
        q11,
        q12,
        q1,
        q2,
        q3,
        q4,
        q5,
        q6,
        q7,
        q8,
        q9,
        q10,
        q11,
        q12,
        q1,
        q2,
        q3,
        q4,
        q5,
        q6,
        q7,
        q8,
        q9,
        q10,
        q11,
        q12,
    ]
).reshape(1, -1)
qq_all = np.repeat(qq, num_refer_points, axis=0)

numtimesteps = 36

# analytical solver
numtemppoints = len(point_x)
T2 = np.zeros([numtemppoints, numtimesteps])

coeff_all = np.zeros([numtemppoints, numtimesteps])

for currstep in range(0, numtimesteps):
    Temp_po_to_referencepo = np.zeros([numtemppoints, numbhe])
    po_dist_to_referencepo = np.zeros([numtemppoints, numbhe])
    localcoeff_all = np.zeros([numtemppoints, 1])
    localcoeff = np.zeros([numtemppoints, numbhe])
    localcoeff1 = np.zeros([numtemppoints, numbhe])
    for i in range(0, numbhe):
        if time_trans * (currstep + 1) - time_trans * 0 > 0:
            for j in range(0, numtemppoints):
                po_dist_to_referencepo[j, i] = (
                    abs(po_x[i] - point_x[j]) ** 2 + abs(po_y[i] - point_y[j]) ** 2
                )
                exp = po_dist_to_referencepo[j, i] / (
                    4 * alpha_sp * time_trans * (currstep + 1)
                )
                n = sp.exp1(exp)
                localcoeff[j, i] = 1 / (4 * math.pi * lamda_sp) * n
        if time_trans * (currstep + 1) - time_trans * 1 > 0:
            for j in range(0, numtemppoints):
                po_dist_to_referencepo[j, i] = (
                    abs(po_x[i] - point_x[j]) ** 2 + abs(po_y[i] - point_y[j]) ** 2
                )
                exp1 = po_dist_to_referencepo[j, i] / (
                    4 * alpha_sp * time_trans * currstep
                )
                n1 = sp.exp1(exp1)
                localcoeff[j, i] = localcoeff[j, i] - 1 / (4 * math.pi * lamda_sp) * n1

    localcoeff_all = np.sum(localcoeff, axis=1).reshape(-1, 1)
    coeff_all[:, 1:] = coeff_all[:, : numtimesteps - 1]
    coeff_all[:, :1] = localcoeff_all

for currstep in range(0, numtimesteps):
    T2[:, currstep] = (
        np.sum(
            coeff_all[:, numtimesteps - 1 - currstep :] * qq_all[:, : currstep + 1],
            axis=1,
        )
        + T0
    )

T2_ini = np.zeros([num_refer_points, 1]) + T0
T2 = np.concatenate((T2_ini, T2), axis=1)

T2_trans = T2

#%% plotting

png_num = 1
for i in range(png_num):
    plt.figure()
    plt.plot(point_x, T2[:, 4], "b", label="Analytical")
    plt.xlim([0, 100])
    plt.ylim([-10, 20])
    plt.ylabel("Temperature [$^\circ$C]")
    plt.xlabel("x [m]")
    plt.legend(loc="best", fontsize=8)
    plt.title(
        f"Soil temperature distribution on A-A'section after 4 months", fontsize=12
    )
    plt.savefig("pngfile{}.png".format(i), dpi=300, transparent=False)
