// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

namespace NumLib
{
template <class T_X, class T_N>
void ShapeTri3::computeShapeFunction(const T_X& r, T_N& N)
{
    N[0] = 1. - r[0] - r[1];
    N[1] = r[0];
    N[2] = r[1];
}

template <class T_X, class T_N>
void ShapeTri3::computeGradShapeFunction(const T_X& /*r*/, T_N& dN)
{
    // dN/dr
    dN[0] = -1.0;
    dN[1] = 1.0;
    dN[2] = 0.0;
    // dN/ds
    dN[3] = -1.0;
    dN[4] = 0.0;
    dN[5] = 1.0;
}

}  // namespace NumLib
