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
void ShapePrism15::computeShapeFunction(const T_X &x, T_N &N)
{
    const double L1 = x[0];
    const double L2 = x[1];
    const double L0 = 1.0 - L1 - L2;
    const double t = x[2];
    const double tt1 = 1.0 - t * t;

    double v1 = 2.0 * L0 - 1;
    double v2 = 2.0 * L1 - 1;
    double v3 = 2.0 * L2 - 1;
    // Vertex, bottom
    N[0] = 0.5 * L0 * (v1 * (1.0 - t) - tt1);
    N[1] = 0.5 * L1 * (v2 * (1.0 - t) - tt1);
    N[2] = 0.5 * L2 * (v3 * (1.0 - t) - tt1);
    // Vertex, top
    N[3] = 0.5 * L0 * (v1 * (1.0 + t) - tt1);
    N[4] = 0.5 * L1 * (v2 * (1.0 + t) - tt1);
    N[5] = 0.5 * L2 * (v3 * (1.0 + t) - tt1);

    v1 = 2.0 * L0 * L1;
    v2 = 2.0 * L1 * L2;
    v3 = 2.0 * L2 * L0;
    // Middle point, bottom
    N[6] = v1 * (1.0 - t);
    N[7] = v2 * (1.0 - t);
    N[8] = v3 * (1.0 - t);
    // Middle point, top
    N[9] = v1 * (1.0 + t);
    N[10] = v2 * (1.0 + t);
    N[11] = v3 * (1.0 + t);
    // Middle point, center
    N[12] = L0 * tt1;
    N[13] = L1 * tt1;
    N[14] = L2 * tt1;
}

template <class T_X, class T_N>
void ShapePrism15::computeGradShapeFunction(const T_X &x, T_N &dN)
{
    const double L1 = x[0];
    const double L2 = x[1];
    const double L0 = 1.0 - L1 - L2;
    const double t = x[2];
    const double tt1 = 1.0 - t * t;

    //---dN/dL1
    double v1 = (4.0 * L0 - 1);
    double v2 = (4.0 * L1 - 1);
    // Vertex, bottom
    dN[0] = -0.5 * (v1 * (1.0 - t) - tt1);
    dN[1] =  0.5 * (v2 * (1.0 - t) - tt1);
    dN[2] =  0.0;
    // Vertex, top
    dN[3] = -0.5 * (v1 * (1.0 + t) - tt1);
    dN[4] =  0.5 * (v2 * (1.0 + t) - tt1);
    dN[5] =  0.0;
    // Middle point, bottom
    dN[6] =  2.0 * (L0 -L1) * (1.0 - t);
    dN[7] =  2.0 * L2 * (1.0 - t);
    dN[8] = -dN[7];
    // Middle point, top
    dN[9] =  2.0 * (L0 - L1) * (1.0 + t);
    dN[10] = 2.0 * L2 * (1.0 + t);
    dN[11] = -dN[10];
    // Middle point, center
    dN[12] = -tt1;
    dN[13] =  tt1;
    dN[14] =  0.0;

    //---dN/dL2
    v1 = (4.0 * L2 - 1);
    // Vertex, bottom
    dN[15] =  dN[0];
    dN[16] =  0.0;
    dN[17] =  0.5 * (v1 * (1.0 - t) - tt1);
    // Vertex, top
    dN[18] =  dN[3];
    dN[19] =  0.0;
    dN[20] =  0.5 * (v1 * (1.0 + t) - tt1);
    // Middle point, bottom
    dN[21] = -2.0 * L1 * (1.0 - t);
    dN[22] = -dN[21];
    v1 = 2.0 * (L0 -  L2);
    dN[23] =  v1 * (1.0 - t);
    // Middle point, top
    dN[24] = -2.0 * L1 * (1.0 + t);
    dN[25] = -dN[24];
    dN[26] =  v1 * (1.0 + t);
    // Middle point, center
    dN[27] = -tt1;
    dN[28] =  0.0;
    dN[29] =  tt1;

    //---dN/dt
    v1 = 2.0 * L0 - 1;
    v2 = 2.0 * L1 - 1;
    double v3 = 2.0 * L2 - 1;
    // Vertex, bottom
    dN[30] = 0.5 * L0 * (-v1 + 2.0 * t);
    dN[31] = 0.5 * L1 * (-v2 + 2.0 * t);
    dN[32] = 0.5 * L2 * (-v3 + 2.0 * t);
    // Vertex, top
    dN[33] = 0.5 * L0 * (v1 + 2.0 * t);
    dN[34] = 0.5 * L1 * (v2 + 2.0 * t);
    dN[35] = 0.5 * L2 * (v3 + 2.0 * t);
    // Middle point, bottom
    dN[36] = -2.0 * L0 * L1;
    dN[37] = -2.0 * L1 * L2;
    dN[38] = -2.0 * L2 * L0;
    // Middle point, top
    dN[39] = -dN[36];
    dN[40] = -dN[37];
    dN[41] = -dN[38];
    // Middle point, center
    dN[42] = -2.0 * L0 * t;
    dN[43] = -2.0 * L1 * t;
    dN[44] = -2.0 * L2 * t;
}

}

