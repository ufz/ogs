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
void ShapeHex8::computeShapeFunction(const T_X &r, T_N &N)
{
    N[0] = (1.0 - r[0]) * (1.0 - r[1]) * (1.0 - r[2]) * 0.125;
    N[1] = (1.0 + r[0]) * (1.0 - r[1]) * (1.0 - r[2]) * 0.125;
    N[2] = (1.0 + r[0]) * (1.0 + r[1]) * (1.0 - r[2]) * 0.125;
    N[3] = (1.0 - r[0]) * (1.0 + r[1]) * (1.0 - r[2]) * 0.125;
    N[4] = (1.0 - r[0]) * (1.0 - r[1]) * (1.0 + r[2]) * 0.125;
    N[5] = (1.0 + r[0]) * (1.0 - r[1]) * (1.0 + r[2]) * 0.125;
    N[6] = (1.0 + r[0]) * (1.0 + r[1]) * (1.0 + r[2]) * 0.125;
    N[7] = (1.0 - r[0]) * (1.0 + r[1]) * (1.0 + r[2]) * 0.125;
}

template <class T_X, class T_N>
void ShapeHex8::computeGradShapeFunction(const T_X &r, T_N &dN)
{
    // dN/dx
    dN[0] = -(1.0 - r[1]) * (1.0 - r[2]) * 0.125;
    dN[1] = -dN[0];
    dN[2] = +(1.0 + r[1]) * (1.0 - r[2]) * 0.125;
    dN[3] = -dN[2];
    dN[4] = -(1.0 - r[1]) * (1.0 + r[2]) * 0.125;
    dN[5] = -dN[4];
    dN[6] = +(1.0 + r[1]) * (1.0 + r[2]) * 0.125;
    dN[7] = -dN[6];

    // dN/dy
    dN[8]  = -(1.0 - r[0]) * (1.0 - r[2]) * 0.125;
    dN[9]  = -(1.0 + r[0]) * (1.0 - r[2]) * 0.125;
    dN[10] = -dN[9];
    dN[11] = -dN[8];
    dN[12] = -(1.0 - r[0]) * (1.0 + r[2]) * 0.125;
    dN[13] = -(1.0 + r[0]) * (1.0 + r[2]) * 0.125;
    dN[14] = -dN[13];
    dN[15] = -dN[12];

    // dN/dz
    dN[16] = -(1.0 - r[0]) * (1.0 - r[1]) * 0.125;
    dN[17] = -(1.0 + r[0]) * (1.0 - r[1]) * 0.125;
    dN[18] = -(1.0 + r[0]) * (1.0 + r[1]) * 0.125;
    dN[19] = -(1.0 - r[0]) * (1.0 + r[1]) * 0.125;
    dN[20] = -dN[16];
    dN[21] = -dN[17];
    dN[22] = -dN[18];
    dN[23] = -dN[19];
}

}

