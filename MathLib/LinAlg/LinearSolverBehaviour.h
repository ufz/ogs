// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

namespace MathLib
{
enum class LinearSolverBehaviour : int
{
    RECOMPUTE,
    RECOMPUTE_AND_STORE,
    REUSE
};

}  // namespace MathLib
