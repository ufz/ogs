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
void ShapeQuad9::computeGradShapeFunction(const T_X& r, T_N& dNdr)
{
    dNdr[0] = (r[0] + 0.5) * r[1] * (r[1] + 1) / 2;
    dNdr[1] = (r[0] - 0.5) * r[1] * (r[1] + 1) / 2;
    dNdr[2] = (r[0] - 0.5) * r[1] * (r[1] - 1) / 2;
    dNdr[3] = (r[0] + 0.5) * r[1] * (r[1] - 1) / 2;
    dNdr[4] = -r[0] * r[1] * (1 + r[1]);
    dNdr[5] = (1 - r[1] * r[1]) * (r[0] - 0.5);
    dNdr[6] = r[0] * r[1] * (1 - r[1]);
    dNdr[7] = (1 - r[1] * r[1]) * (r[0] + 0.5);
    dNdr[8] = 2 * r[0] * (r[1] * r[1] - 1);

    dNdr[10] = (r[1] + 0.5) * r[0] * (r[0] - 1) / 2;
    dNdr[11] = (r[1] - 0.5) * r[0] * (r[0] - 1) / 2;
    dNdr[12] = (r[1] - 0.5) * r[0] * (r[0] + 1) / 2;
    dNdr[13] = (1 - r[0] * r[0]) * (r[1] + 0.5);
    dNdr[14] = r[0] * r[1] * (1 - r[0]);
    dNdr[15] = (1 - r[0] * r[0]) * (r[1] - 0.5);
    dNdr[16] = -r[0] * r[1] * (1 + r[0]);
    dNdr[17] = 2 * r[1] * (r[0] * r[0] - 1);
    dNdr[9] = (r[1] + 0.5) * r[0] * (r[0] + 1) / 2;
}
}

