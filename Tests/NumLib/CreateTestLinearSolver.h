// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>

#if defined(USE_LIS)
#include "MathLib/LinAlg/EigenLis/LinearSolverOptionsParser.h"
#elif defined(USE_PETSC)
#include "MathLib/LinAlg/PETSc/LinearSolverOptionsParser.h"
#else
#include "MathLib/LinAlg/Eigen/LinearSolverOptionsParser.h"
#endif
#include "MathLib/LinAlg/GlobalLinearSolverType.h"

//! Creates the default linear solver used by the ODE system tests.
inline std::unique_ptr<GlobalLinearSolver> createLinearSolver()
{
#if defined(USE_PETSC)
    std::string const petsc_options =
        "-ksp_type bcgs -pc_type sor -ksp_rtol 1e-24 -ksp_max_it 100 "
        "-ksp_initial_guess_nonzero false";
    return std::make_unique<GlobalLinearSolver>("", petsc_options);
#else
    auto const solver_options =
        MathLib::LinearSolverOptionsParser<GlobalLinearSolver>{}
            .parseNameAndOptions("", nullptr);
    return std::make_unique<GlobalLinearSolver>(std::get<0>(solver_options),
                                                std::get<1>(solver_options));
#endif
}
