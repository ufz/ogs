/**
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
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
    N[0] = (1.0 - r[0]) * (1.0 - r[1]) * (1.0 - r[2]);
    N[1] = (1.0 + r[0]) * (1.0 - r[1]) * (1.0 - r[2]);
    N[2] = (1.0 + r[0]) * (1.0 + r[1]) * (1.0 - r[2]);
    N[3] = (1.0 - r[0]) * (1.0 + r[1]) * (1.0 - r[2]);
    N[4] = (1.0 - r[0]) * (1.0 - r[1]) * (1.0 + r[2]);
    N[5] = (1.0 + r[0]) * (1.0 - r[1]) * (1.0 + r[2]);
    N[6] = (1.0 + r[0]) * (1.0 + r[1]) * (1.0 + r[2]);
    N[7] = (1.0 - r[0]) * (1.0 + r[1]) * (1.0 + r[2]);
    for (unsigned i = 0; i < 8; i++)
        N[i] *= 0.125;
}

template <class T_X, class T_N>
void ShapeHex8::computeGradShapeFunction(const T_X &r, T_N &dN)
{
    // dN/dx
    dN[0] = -(1.0 - r[1]) * (1.0 - r[2]);
    dN[1] = +(1.0 - r[1]) * (1.0 - r[2]);
    dN[2] = +(1.0 + r[1]) * (1.0 - r[2]);
    dN[3] = -(1.0 + r[1]) * (1.0 - r[2]);
    dN[4] = -(1.0 - r[1]) * (1.0 + r[2]);
    dN[5] = +(1.0 - r[1]) * (1.0 + r[2]);
    dN[6] = +(1.0 + r[1]) * (1.0 + r[2]);
    dN[7] = -(1.0 + r[1]) * (1.0 + r[2]);

    // dN/dy
    dN[8]  = -(1.0 - r[0]) * (1.0 - r[2]);
    dN[9]  = -(1.0 + r[0]) * (1.0 - r[2]);
    dN[10] = +(1.0 + r[0]) * (1.0 - r[2]);
    dN[11] = +(1.0 - r[0]) * (1.0 - r[2]);
    dN[12] = -(1.0 - r[0]) * (1.0 + r[2]);
    dN[13] = -(1.0 + r[0]) * (1.0 + r[2]);
    dN[14] = +(1.0 + r[0]) * (1.0 + r[2]);
    dN[15] = +(1.0 - r[0]) * (1.0 + r[2]);

    // dN/dz
    dN[16] = -(1.0 - r[0]) * (1.0 - r[1]);
    dN[17] = -(1.0 + r[0]) * (1.0 - r[1]);
    dN[18] = -(1.0 + r[0]) * (1.0 + r[1]);
    dN[19] = -(1.0 - r[0]) * (1.0 + r[1]);
    dN[20] = +(1.0 - r[0]) * (1.0 - r[1]);
    dN[21] = +(1.0 + r[0]) * (1.0 - r[1]);
    dN[22] = +(1.0 + r[0]) * (1.0 + r[1]);
    dN[23] = +(1.0 - r[0]) * (1.0 + r[1]);

    for (unsigned i = 0; i < 24; i++)
        dN[i] *= 0.125;
}

}

