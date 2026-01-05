// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

namespace MathLib
{
template <typename MAT_T>
bool finalizeMatrixAssembly(MAT_T& /*unused*/)
{
    return true;
}

}  // namespace MathLib
