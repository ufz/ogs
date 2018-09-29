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
void ShapeLine3::computeShapeFunction(const T_X &r, T_N &N)
{
    N[0] = 0.5 * r[0] * (r[0] - 1.0);
    N[1] = 0.5 * r[0] * (r[0] + 1.0);
    N[2] = 1.0 - r[0] * r[0];
}

template <class T_X, class T_N>
void ShapeLine3::computeGradShapeFunction(const T_X &r, T_N &dN)
{
    dN[0] = r[0] - 0.5;
    dN[1] = r[0] + 0.5;
    dN[2] = -2.0 * r[0];
}

}

