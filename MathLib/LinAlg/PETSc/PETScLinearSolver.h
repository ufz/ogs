// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <petscksp.h>

#include <algorithm>
#include <array>
#include <string>

#include "MathLib/LinAlg/LinearSolverBehaviour.h"
#include "PETScMatrix.h"
#include "PETScVector.h"

namespace MathLib
{
/*! A class of linear solver based on PETSc routines.

 All command-line options that are not recognized by OGS are passed on to
 PETSc, i.e., they potentially affect the linear solver.
 The linear solver options in the project file take precedence over the
 command-line options, because the former are processed at a later stage.
*/
class PETScLinearSolver
{
public:
    /*!
        Constructor
        \param prefix  Name used to distinguish the options in the command
                       line for this solver. It can be the name of the PDE
                       that owns an instance of this class.
        \param petsc_options PETSc options string which is passed to PETSc lib
        and inserted in the PETSc option database (see
        https://petsc.org/release/manualpages/Sys/PetscOptionsInsertString/).
    */
    PETScLinearSolver(std::string const& prefix,
                      std::string const& petsc_options);

    ~PETScLinearSolver() { KSPDestroy(&solver_); }
    // TODO check if some args in LinearSolver interface can be made const&.
    bool solve(PETScMatrix& A, PETScVector& b, PETScVector& x);

    /// Get number of iterations.
    PetscInt getNumberOfIterations() const
    {
        PetscInt its = 0;
        KSPGetIterationNumber(solver_, &its);
        return its;
    }

    /// Get elapsed wall clock time.
    double getElapsedTime() const { return elapsed_ctime_; }

    /// Get, if the solver can handle rectangular equation systems
    bool canSolveRectangular() const { return can_solve_rectangular_; }

    /// Provided for compatibility with Eigen linear solvers.
    bool willCompute(
        MathLib::LinearSolverBehaviour const /*linear_solver_behaviour*/) const
    {
        return true;
    }

private:
    bool canSolverHandleRectangular(std::string_view ksp_type)
    {
        // List of KSP types that can handle rectangular matrices
        constexpr auto rectangular_solvers = std::to_array<std::string_view>(
            {KSPGMRES, KSPFGMRES, KSPBCGS, KSPBCGSL, KSPTFQMR, KSPCGS});

        return std::ranges::any_of(rectangular_solvers,
                                   [&](const auto& solver)
                                   { return solver == ksp_type; });
    }

    KSP solver_;  ///< Solver type.
    PC pc_;       ///< Preconditioner type.

    bool can_solve_rectangular_ = false;

    double elapsed_ctime_ = 0.0;  ///< Clock time
};

}  // namespace MathLib
