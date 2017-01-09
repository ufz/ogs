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
void ShapePrism6::computeShapeFunction(const T_X &x, T_N &N)
{
    const double L1 = x[0];
    const double L2 = x[1];
    const double t = x[2];
    N[0] = 0.5 * (1.0 - L1 - L2) * (1.0 - t);
    N[1] = 0.5 * L1 * (1.0 - t);
    N[2] = 0.5 * L2 * (1.0 - t);
    N[3] = 0.5 * (1.0 - L1 - L2) * (1.0 + t);
    N[4] = 0.5 * L1 * (1.0 + t);
    N[5] = 0.5 * L2 * (1.0 + t);
}

template <class T_X, class T_N>
void ShapePrism6::computeGradShapeFunction(const T_X &x, T_N &dN)
{
    const double L1 = x[0];
    const double L2 = x[1];
    const double t = x[2];
    //  dN/dL1
    dN[0] = -0.5 * (1.0 - t);
    dN[1] = -dN[0];
    dN[2] =  0.0;
    dN[3] = -0.5 * (1.0 + t);
    dN[4] = -dN[3];
    dN[5] =  0.0;
    //  dN/dL2
    dN[6] =  dN[0];
    dN[7] =  0.0;
    dN[8] = -dN[0];
    dN[9] =  dN[3];
    dN[10] = 0.0;
    dN[11] = -dN[3];
    //  dN/dt
    dN[12] = -0.5 * (1.0 - L1 - L2);
    dN[13] = -0.5 * L1;
    dN[14] = -0.5 * L2;
    dN[15] = -dN[12];
    dN[16] = -dN[13];
    dN[17] = -dN[14];
}

}

