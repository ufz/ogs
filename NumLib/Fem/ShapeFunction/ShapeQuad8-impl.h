/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

namespace NumLib
{

template <class T_X, class T_N>
void ShapeQuad8::computeShapeFunction(const T_X &r, T_N &N)
{
    N[0] = 0.25 * (1.0 + r[0]) * (1.0 + r[1]) * (-1.0 + r[0] + r[1]);
    N[1] = -0.25 * (1.0 - r[0]) * (1.0 + r[1]) * (1.0 + r[0] - r[1]);
    N[2] = -0.25 * (1.0 - r[0]) * (1.0 - r[1]) * (1.0 + r[0] + r[1]);
    N[3] = 0.25 * (1.0 + r[0]) * (1.0 - r[1]) * (-1.0 + r[0] - r[1]);
    //
    N[4] = 0.5 * (1.0 - r[0] * r[0]) * (1.0 + r[1]);
    N[5] = 0.5 * (1.0 - r[1] * r[1]) * (1.0 - r[0]);
    N[6] = 0.5 * (1.0 - r[0] * r[0]) * (1.0 - r[1]);
    N[7] = 0.5 * (1.0 - r[1] * r[1]) * (1.0 + r[0]);
}

template <class T_X, class T_N>
void ShapeQuad8::computeGradShapeFunction(const T_X &rs, T_N &dNdr)
{
    const double r = rs[0];
    const double s = rs[1];

    // dN/dr
    dNdr[0] = (1 + s) * (2 * r + s) * 0.25;
    dNdr[1] = (1 + s) * (2 * r - s) * 0.25;
    dNdr[2] = (1 - s) * (2 * r + s) * 0.25;
    dNdr[3] = (1 - s) * (2 * r - s) * 0.25;

    dNdr[4] = -r * (1 + s);
    dNdr[5] = -(1 - s * s) * 0.5;
    dNdr[6] = -r * (1 - s);
    dNdr[7] = (1 - s * s) * 0.5;

    // dN/ds
    dNdr[8] = (1 + r) * (r + 2 * s) * 0.25;
    dNdr[9] = -(1 - r) * (r - 2 * s) * 0.25;
    dNdr[10] = (1 - r) * (r + 2 * s) * 0.25;
    dNdr[11] = -(1 + r) * (r - 2 * s) * 0.25;

    dNdr[12] = (1 - r * r) * 0.5;
    dNdr[13] = -(1 - r) * s;
    dNdr[14] = -(1 - r * r) * 0.5;
    dNdr[15] = -(1 + r) * s;
}

}

