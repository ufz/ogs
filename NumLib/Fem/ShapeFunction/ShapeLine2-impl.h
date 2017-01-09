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
void ShapeLine2::computeShapeFunction(const T_X &r, T_N &N)
{
    N[0] = (1.0 - r[0]) * 0.5;
    N[1] = (1.0 + r[0]) * 0.5;
}

template <class T_X, class T_N>
void ShapeLine2::computeGradShapeFunction(const T_X &/*r*/, T_N &dN)
{
    dN[0] = -0.5;
    dN[1] = 0.5;
}

}

