/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once
#ifdef USE_PETSC

#include <petscsnes.h>

#include "ConvergenceCriterion.h"
#include "NonlinearSolver.h"
#include "TimeDiscretizedODESystem.h"

namespace NumLib
{
//! \addtogroup ODESolver
//! @{
class PETScNonlinearSolver final : public NonlinearSolverBase
{
public:
    //! Type of the nonlinear equation system to be solved.
    using System = NonlinearSystem<NonlinearSolverTag::Newton>;

    /*! Constructs a new instance.
     *
     * \param linear_solver the linear solver used by this nonlinear solver.
     * \param maxiter the maximum number of iterations used to solve the
     *                equation.
     * \param prefix used to set the SNES options.
     */
    PETScNonlinearSolver(GlobalLinearSolver& linear_solver,
                         int maxiter,
                         std::string prefix);

    //! Set the nonlinear equation system that will be solved.
    void setEquationSystem(System& eq, ConvergenceCriterion& conv_crit);

    void compensateNonEquilibriumInitialResiduum(bool const value);

    void calculateNonEquilibriumInitialResiduum(
        std::vector<GlobalVector*> const& x,
        std::vector<GlobalVector*> const& x_prev,
        int const process_id) override;

    NonlinearSolverStatus solve(
        std::vector<GlobalVector*>& x,
        std::vector<GlobalVector*> const& x_prev,
        std::function<void(int, std::vector<GlobalVector*> const&)> const&
            postIterationCallback,
        int const process_id) override;

private:
    SNES _snes_solver;

    System* _equation_system = nullptr;
    ConvergenceCriterion* _convergence_criterion = nullptr;

    std::size_t _residual_id = 0u;  //!< ID of the residual vector.
    std::size_t _jacobian_id = 0u;  //!< ID of the Jacobian matrix.

    std::size_t _petsc_x_id = 0u;
    std::size_t _petsc_jacobian_id = 0u;
    std::size_t _petsc_residual_id = 0u;

    // clang-format off
    /// \copydoc NumLib::NonlinearSolver<NonlinearSolverTag::Newton>::_compensate_non_equilibrium_initial_residuum
    bool _compensate_non_equilibrium_initial_residuum = false;
    // clang-format on
};

//! @}

}  // namespace NumLib
#endif  // USE_PETSC
