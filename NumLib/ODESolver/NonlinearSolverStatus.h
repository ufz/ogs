// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

namespace NumLib
{
/// Status of the non-linear solver.
struct NonlinearSolverStatus
{
    bool error_norms_met = false;
    int number_iterations = -1;
};
}  // namespace NumLib
