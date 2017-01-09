/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

namespace NumLib
{

template <class T_X, class T_N>
void ShapeQuad9::computeShapeFunction(const T_X &r, T_N &N)
{
    N[8] = (1.0 - r[0] * r[0]) * ( 1.0 - r[1] * r[1]);
    N[7] = 0.5 * (1.0 - r[1] * r[1]) * (1.0 + r[0]) - 0.5 * N[8];
    N[6] = 0.5 * (1.0 - r[0] * r[0]) * (1.0 - r[1]) - 0.5 * N[8];
    N[5] = 0.5 * (1.0 - r[1] * r[1]) * (1.0 - r[0]) - 0.5 * N[8];
    N[4] = 0.5 * (1.0 - r[0] * r[0]) * (1.0 + r[1]) - 0.5 * N[8];
    N[3] = 0.25 * (1.0 + r[0]) * (1.0 - r[1]) - 0.5 * N[6] - 0.5 * N[7] - 0.25 * N[8];
    N[2] = 0.25 * (1.0 - r[0]) * (1.0 - r[1]) - 0.5 * N[5] - 0.5 * N[6] - 0.25 * N[8];
    N[1] = 0.25 * (1.0 - r[0]) * (1.0 + r[1]) - 0.5 * N[4] - 0.5 * N[5] - 0.25 * N[8];
    N[0] = 0.25 * (1.0 + r[0]) * (1.0 + r[1]) - 0.5 * N[4] - 0.5 * N[7] - 0.25 * N[8];
}

template <class T_X, class T_N>
void ShapeQuad9::computeGradShapeFunction(const T_X &r, T_N &dNdr)
{
    dNdr[8] = -2.0 * r[0] * (1.0 - r[1] * r[1]);
    dNdr[7] = +0.5 * (1.0 - r[1] * r[1]) - 0.5 * dNdr[8];
    dNdr[6] = -1.0 * r[0] * (1.0 - r[1]) - 0.5 * dNdr[8];
    dNdr[5] = -dNdr[7];
    dNdr[4] = -1.0 * r[0] * (1.0 + r[1]) - 0.5 * dNdr[8];
    dNdr[3] = +0.25 * (1 - r[1]) - 0.5 * dNdr[6] - 0.5 * dNdr[7] - 0.25 * dNdr[8];
    dNdr[2] = -0.25 * (1 - r[1]) - 0.5 * dNdr[5] - 0.5 * dNdr[6] - 0.25 * dNdr[8];
    dNdr[1] = -0.25 * (1 + r[1]) - 0.5 * dNdr[4] - 0.5 * dNdr[5] - 0.25 * dNdr[8];
    dNdr[0] = +0.25 * (1 + r[1]) - 0.5 * dNdr[4] - 0.5 * dNdr[7] - 0.25 * dNdr[8];

    dNdr[17] = -2.0 * r[1] * (1.0 - r[0] * r[0]);
    dNdr[16] = -1.0 * r[1] * (1.0 + r[0]) - 0.5 * dNdr[17];
    dNdr[15] = -0.5 * (1.0 - r[0] * r[0]) - 0.5 * dNdr[17];
    dNdr[14] = -1.0 * r[1] * (1.0 - r[0]) - 0.5 * dNdr[17];
    dNdr[13] = -dNdr[15];
    dNdr[12] = -0.25 * (1 + r[0]) - 0.5 * dNdr[15] - 0.5 * dNdr[16] - 0.25 * dNdr[17];
    dNdr[11] = -0.25 * (1 - r[0]) - 0.5 * dNdr[14] - 0.5 * dNdr[15] - 0.25 * dNdr[17];
    dNdr[10] = +0.25 * (1 - r[0]) - 0.5 * dNdr[13] - 0.5 * dNdr[14] - 0.25 * dNdr[17];
    dNdr[9] = +0.25 * (1 + r[0]) - 0.5 * dNdr[13] - 0.5 * dNdr[16] - 0.25 * dNdr[17];
}

}

