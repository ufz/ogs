/**
 * \author Norihiro Watanabe
 * \date   2013-08-13
 *
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
void ShapeQuad4::computeShapeFunction(const T_X &r, T_N &N)
{
    N[0] = (1.0 + r[0]) * (1.0 + r[1]);
    N[1] = (1.0 - r[0]) * (1.0 + r[1]);
    N[2] = (1.0 - r[0]) * (1.0 - r[1]);
    N[3] = (1.0 + r[0]) * (1.0 - r[1]);
    for (int i = 0; i < 4; i++)
        N[i] *= 0.25;
}

template <class T_X, class T_N>
void ShapeQuad4::computeGradShapeFunction(const T_X &r, T_N &dN)
{
    dN[0] = +(1.0 + r[1]);
    dN[1] = -(1.0 + r[1]);
    dN[2] = -(1.0 - r[1]);
    dN[3] = +(1.0 - r[1]);
    dN[4] = +(1.0 + r[0]);
    dN[5] = +(1.0 - r[0]);
    dN[6] = -(1.0 - r[0]);
    dN[7] = -(1.0 + r[0]);
    for (int i = 0; i < 8; i++)
        dN[i] *= 0.25;
}

}

