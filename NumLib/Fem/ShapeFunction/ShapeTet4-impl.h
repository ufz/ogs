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
void ShapeTet4::computeShapeFunction(const T_X &r, T_N &N)
{
    N[0] = 1. - r[0] - r[1] - r[2];
    N[1] = r[0];
    N[2] = r[1];
    N[3] = r[2];
}

template <class T_X, class T_N>
void ShapeTet4::computeGradShapeFunction(const T_X &/*r*/, T_N &dNdr)
{
    //dr
    dNdr[0] = -1.0;
    dNdr[1] = 1.0;
    dNdr[2] = 0.0;
    dNdr[3] = 0.0;

    //ds
    dNdr[4] = -1.0;
    dNdr[5] = 0.0;
    dNdr[6] = 1.0;
    dNdr[7] = 0.0;

    //dt
    dNdr[8] = -1.0;
    dNdr[9] = 0.0;
    dNdr[10] = 0.0;
    dNdr[11] = 1.0;
}

}

