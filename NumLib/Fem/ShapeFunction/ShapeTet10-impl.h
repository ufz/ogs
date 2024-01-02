/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

namespace NumLib
{
template <class T_X, class T_N>
void ShapeTet10::computeShapeFunction(const T_X& r, T_N& N)
{
    N[0] = 2. * (1 - r[0] - r[1] - r[2]) * (0.5 - r[0] - r[1] - r[2]);
    N[1] = r[0] * (2. * r[0] - 1);
    N[2] = r[1] * (2. * r[1] - 1);
    N[3] = r[2] * (2. * r[2] - 1);
    N[4] = 4.0 * r[0] * (1.0 - r[0] - r[1] - r[2]);
    N[5] = 4.0 * r[0] * r[1];
    N[6] = 4.0 * r[1] * (1.0 - r[0] - r[1] - r[2]);
    N[7] = 4.0 * r[2] * (1.0 - r[0] - r[1] - r[2]);
    N[8] = 4.0 * r[0] * r[2];
    N[9] = 4.0 * r[1] * r[2];
}

template <class T_X, class T_N>
void ShapeTet10::computeGradShapeFunction(const T_X& r, T_N& dN)
{
    dN[0] = 4.0 * (r[0] + r[1] + r[2]) - 3.0;
    dN[1] = 4. * r[0] - 1.;
    dN[2] = 0.0;
    dN[3] = 0.0;
    dN[4] = 4.0 * (1.0 - 2.0 * r[0] - r[1] - r[2]);
    dN[5] = 4.0 * r[1];
    dN[6] = -4.0 * r[1];
    dN[7] = -4.0 * r[2];
    dN[8] = 4.0 * r[2];
    dN[9] = 0.0;

    dN[10] = 4. * (r[0] + r[1] + r[2]) - 3.;
    dN[11] = 0.0;
    dN[12] = 4. * r[1] - 1.;
    dN[13] = 0.;
    dN[14] = -4.0 * r[0];
    dN[15] = 4.0 * r[0];
    dN[16] = 4.0 * (1.0 - r[0] - 2.0 * r[1] - r[2]);
    dN[17] = -4.0 * r[2];
    dN[18] = 0.0;
    dN[19] = 4.0 * r[2];

    dN[20] = 4. * (r[0] + r[1] + r[2]) - 3.;
    dN[21] = 0.;
    dN[22] = 0.;
    dN[23] = 4. * r[2] - 1.;
    dN[24] = -4.0 * r[0];
    dN[25] = 0.0;
    dN[26] = -4.0 * r[1];
    dN[27] = 4.0 * (1.0 - r[0] - r[1] - 2.0 * r[2]);
    dN[28] = 4.0 * r[0];
    dN[29] = 4.0 * r[1];
}

}  // namespace NumLib
