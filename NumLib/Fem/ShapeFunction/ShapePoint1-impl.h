// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

namespace NumLib
{
template <class T_X, class T_N>
void ShapePoint1::computeShapeFunction(const T_X& /*r*/, T_N& N)
{
    N[0] = 1;
}

}  // namespace NumLib
