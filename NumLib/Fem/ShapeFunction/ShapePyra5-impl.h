// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

namespace NumLib
{
template <class T_X, class T_N>
void ShapePyra5::computeShapeFunction(const T_X& x, T_N& N)
{
    const double r = x[0];
    const double s = x[1];
    const double t = x[2];

    N[0] = 0.125 * (1 - r) * (1 - s) * (1 - t);
    N[1] = 0.125 * (1 + r) * (1 - s) * (1 - t);
    N[2] = 0.125 * (1 + r) * (1 + s) * (1 - t);
    N[3] = 0.125 * (1 - r) * (1 + s) * (1 - t);
    N[4] = 0.5 * (1 + t);
}

template <class T_X, class T_N>
void ShapePyra5::computeGradShapeFunction(const T_X& x, T_N& dN)
{
    const double r = x[0];
    const double s = x[1];
    const double t = x[2];
    //  dN/dL1
    dN[0] = -0.125 * (1.0 - s) * (1.0 - t);
    dN[1] = -dN[0];
    dN[2] = 0.125 * (1.0 + s) * (1.0 - t);
    dN[3] = -dN[2];
    dN[4] = 0.0;
    //  dN/dL2
    dN[5] = -0.125 * (1.0 - r) * (1.0 - t);
    dN[6] = -0.125 * (1.0 + r) * (1.0 - t);
    dN[7] = -dN[6];
    dN[8] = -dN[5];
    dN[9] = 0.0;
    //  dN/dt
    dN[10] = -0.125 * (1.0 - r) * (1.0 - s);
    dN[11] = -0.125 * (1.0 + r) * (1.0 - s);
    dN[12] = -0.125 * (1.0 + r) * (1.0 + s);
    dN[13] = -0.125 * (1.0 - r) * (1.0 + s);
    dN[14] = 0.5;
}

}  // namespace NumLib
