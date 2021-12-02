/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ProcessData.h"

#include "NumLib/ODESolver/PETScNonlinearSolver.h"

namespace ProcessLib
{
void setEquationSystem(NumLib::NonlinearSolverBase& nonlinear_solver,
                       NumLib::EquationSystem& eq_sys,
                       NumLib::ConvergenceCriterion& conv_crit,
                       NumLib::NonlinearSolverTag nl_tag)
{
    using Tag = NumLib::NonlinearSolverTag;
    switch (nl_tag)
    {
        case Tag::Picard:
        {
            using EqSys = NumLib::NonlinearSystem<Tag::Picard>;
            auto& eq_sys_ = static_cast<EqSys&>(eq_sys);
            if (auto* nl_solver =
                    dynamic_cast<NumLib::NonlinearSolver<Tag::Picard>*>(
                        &nonlinear_solver);
                nl_solver != nullptr)
            {
                nl_solver->setEquationSystem(eq_sys_, conv_crit);
            }
            else
            {
                OGS_FATAL(
                    "Could not cast nonlinear solver to Picard type solver.");
            }
            break;
        }
        case Tag::Newton:
        {
            using EqSys = NumLib::NonlinearSystem<Tag::Newton>;
            auto& eq_sys_ = static_cast<EqSys&>(eq_sys);

            if (auto* nl_solver =
                    dynamic_cast<NumLib::NonlinearSolver<Tag::Newton>*>(
                        &nonlinear_solver);
                nl_solver != nullptr)
            {
                nl_solver->setEquationSystem(eq_sys_, conv_crit);
            }
#ifdef USE_PETSC
            else if (auto* nl_solver =
                         dynamic_cast<NumLib::PETScNonlinearSolver*>(
                             &nonlinear_solver);
                     nl_solver != nullptr)
            {
                nl_solver->setEquationSystem(eq_sys_, conv_crit);
            }
#endif  // USE_PETSC
            else
            {
                OGS_FATAL(
                    "Could not cast nonlinear solver to Newton type solver.");
            }
            break;
        }
    }
}

}  // namespace ProcessLib
