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
void ShapeTri6::computeShapeFunction(const T_X &r, T_N &N)
{
    N[0] = 2. * (1. - r[0] - r[1]) * (0.5 - r[0] - r[1]);
    N[1] = r[0] * (2. * r[0] - 1.);
    N[2] = r[1] * (2. * r[1] - 1.);
    N[3] = 4. * r[0] * (1. - r[0] - r[1]);
    N[4] = 4. * r[0] * r[1];
    N[5] = 4. * r[1] * (1. - r[0] - r[1]);
}

template <class T_X, class T_N>
void ShapeTri6::computeGradShapeFunction(const T_X &r, T_N &dNdr)
{
    dNdr[0] = 4. * (r[0] + r[1]) - 3.;     // dN1/dL1
    dNdr[6] = dNdr[0];                     // dN1/dL2

    dNdr[1] = 4. * r[0] - 1.;              // dN2/dL1
    dNdr[7] = 0.;                          // dN2/dL2

    dNdr[2] = 0.;                          // dN3/dL1
    dNdr[8] = 4. * r[1] - 1.;              // dN3/dL2

    dNdr[3] =  4. * (1 - 2. * r[0] - r[1]); // dN4/dL1
    dNdr[9] = -4. * r[0];                  // dN4/dL2

    dNdr[4] = 4. * r[1];                   // dN5/dL1
    dNdr[10] = -dNdr[9];                   // dN5/dL2

    dNdr[5] = -dNdr[4];                  // dN6/dL1
    dNdr[11] = 4. * (1 - r[0] - 2. * r[1]); // dN6/dL2
}

}

