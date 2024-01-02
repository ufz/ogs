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
void ShapeQuad9::computeShapeFunction(const T_X& r, T_N& N)
{
    N[0] = r[0] * (r[0] + 1) * r[1] * (r[1] + 1) / 4;
    N[1] = r[0] * (r[0] - 1) * r[1] * (r[1] + 1) / 4;
    N[2] = r[0] * (r[0] - 1) * r[1] * (r[1] - 1) / 4;
    N[3] = r[0] * (r[0] + 1) * r[1] * (r[1] - 1) / 4;
    N[4] = r[1] * (r[1] + 1) * (1 - r[0] * r[0]) / 2;
    N[5] = r[0] * (r[0] - 1) * (1 - r[1] * r[1]) / 2;
    N[6] = r[1] * (r[1] - 1) * (1 - r[0] * r[0]) / 2;
    N[7] = r[0] * (r[0] + 1) * (1 - r[1] * r[1]) / 2;
    N[8] = (1 - r[0] * r[0]) * (1 - r[1] * r[1]);
}

template <class T_X, class T_N>
void ShapeQuad9::computeGradShapeFunction(const T_X& r, T_N& dN)
{
    dN[0] = (r[0] + 0.5) * r[1] * (r[1] + 1) / 2;
    dN[1] = (r[0] - 0.5) * r[1] * (r[1] + 1) / 2;
    dN[2] = (r[0] - 0.5) * r[1] * (r[1] - 1) / 2;
    dN[3] = (r[0] + 0.5) * r[1] * (r[1] - 1) / 2;
    dN[4] = -r[0] * r[1] * (1 + r[1]);
    dN[5] = (1 - r[1] * r[1]) * (r[0] - 0.5);
    dN[6] = r[0] * r[1] * (1 - r[1]);
    dN[7] = (1 - r[1] * r[1]) * (r[0] + 0.5);
    dN[8] = 2 * r[0] * (r[1] * r[1] - 1);

    dN[10] = (r[1] + 0.5) * r[0] * (r[0] - 1) / 2;
    dN[11] = (r[1] - 0.5) * r[0] * (r[0] - 1) / 2;
    dN[12] = (r[1] - 0.5) * r[0] * (r[0] + 1) / 2;
    dN[13] = (1 - r[0] * r[0]) * (r[1] + 0.5);
    dN[14] = r[0] * r[1] * (1 - r[0]);
    dN[15] = (1 - r[0] * r[0]) * (r[1] - 0.5);
    dN[16] = -r[0] * r[1] * (1 + r[0]);
    dN[17] = 2 * r[1] * (r[0] * r[0] - 1);
    dN[9] = (r[1] + 0.5) * r[0] * (r[0] + 1) / 2;
}
}  // namespace NumLib
