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
void ShapeTet10::computeShapeFunction(const T_X &r, T_N &N)
{
    N[0] = 2. * (1 - r[0] - r[1] - r[2]) * (0.5 - r[0] - r[1] - r[2]);
    N[1] = r[0] * (2. * r[0] - 1);
    N[2] = r[1] * (2. * r[1] - 1);
    N[3] = r[2] * (2. * r[2] - 1);
    N[4] = 4.0 * r[0] * (1.0 - r[0] - r[1] - r[2]);
    N[5] = 4.0 * r[0] * r[1];
    N[6] = 4.0 * r[1] * (1.0 - r[0] - r[1] - r[2]);
    N[7] = 4.0 * r[0] * r[2];
    N[8] = 4.0 * r[1] * r[2];
    N[9] = 4.0 * r[2] * (1.0 - r[0] - r[1] - r[2]);
}

template <class T_X, class T_N>
void ShapeTet10::computeGradShapeFunction(const T_X &r, T_N &dNdr)
{
    dNdr[0] = 4.0 * (r[0] + r[1] + r[2]) - 3.0;
    dNdr[1] = 4. * r[0] - 1.;
    dNdr[2] = 0.0;
    dNdr[3] = 0.0;
    dNdr[4] = 4.0 * (1.0 - 2.0 * r[0] - r[1] - r[2]);
    dNdr[5] = 4.0 * r[1];
    dNdr[6] = -4.0 * r[1];
    dNdr[7] = 4.0 * r[2];
    dNdr[8] = 0.0;
    dNdr[9] = -4.0 * r[2];

    dNdr[10] =  4. * (r[0] + r[1] + r[2]) - 3.;
    dNdr[11] = 0.0;
    dNdr[12] = 4. * r[1] - 1.;
    dNdr[13] = 0.;
    dNdr[14] = -4.0 * r[0];
    dNdr[15] = 4.0 * r[0];
    dNdr[16] = 4.0 * (1.0 - r[0] - 2.0 * r[1] - r[2]);
    dNdr[17] = 0.0;
    dNdr[18] = 4.0 * r[2];
    dNdr[19] = -4.0 * r[2];

    dNdr[20] = 4. * (r[0] + r[1] + r[2]) - 3.;
    dNdr[21] = 0.;
    dNdr[22] = 0.;
    dNdr[23] = 4. * r[2] - 1.;
    dNdr[24] = -4.0 * r[0];
    dNdr[25] = 0.0;
    dNdr[26] = -4.0 * r[1];
    dNdr[27] = 4.0 * r[0];
    dNdr[28] = 4.0 * r[1];
    dNdr[29] = 4.0 * (1.0 - r[0] - r[1] - 2.0 * r[2]);
}

}

