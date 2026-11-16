// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "ProcessData.h"

#include "NumLib/ODESolver/PETScNonlinearSolver.h"
#include "NumLib/ODESolver/TimeDiscretizedODESystem.h"

namespace ProcessLib
{
template <NumLib::ODESystemTag ODETag>
static void setTimeDiscretizedODESystem(
    ProcessData& process_data,
    NumLib::ODESystem<ODETag, NumLib::NonlinearSolverTag::Picard>& ode_sys,
    ProcessData const* const storage_owner)
{
    auto const make_tdisc_ode_sys = [&]<NumLib::NonlinearSolverTag NLTag>(
                                        NumLib::ODESystem<ODETag, NLTag>& ode)
    {
        using TDiscODESys = NumLib::TimeDiscretizedODESystem<ODETag, NLTag>;
        if (storage_owner == nullptr)
        {
            return std::make_unique<TDiscODESys>(process_data.process_id, ode,
                                                 *process_data.time_disc);
        }
        INFO(
            "The equation system of process {:d} shares the global matrices "
            "and vectors with that of process {:d}.",
            process_data.process_id, storage_owner->process_id);
        return std::make_unique<TDiscODESys>(
            process_data.process_id, ode, *process_data.time_disc,
            dynamic_cast<TDiscODESys const&>(*storage_owner->tdisc_ode_sys));
    };

    using Tag = NumLib::NonlinearSolverTag;
    // A concrete Picard solver
    using NonlinearSolverPicard = NumLib::NonlinearSolver<Tag::Picard>;
    // A concrete Newton solver
    using NonlinearSolverNewton = NumLib::NonlinearSolver<Tag::Newton>;

    if (dynamic_cast<NonlinearSolverPicard*>(&process_data.nonlinear_solver))
    {
        // The Picard solver can also work with a Newton-ready ODE,
        // because the Newton ODESystem derives from the Picard ODESystem.
        // So no further checks are needed here.

        process_data.tdisc_ode_sys = make_tdisc_ode_sys(ode_sys);
    }
    // TODO (naumov) Provide a function to nonlinear_solver to distinguish the
    // types. Could be handy, because a nonlinear solver could handle both types
    // like PETScSNES.
    else if ((dynamic_cast<NonlinearSolverNewton*>(
                  &process_data.nonlinear_solver) != nullptr)
#ifdef USE_PETSC
             || (dynamic_cast<NumLib::PETScNonlinearSolver*>(
                     &process_data.nonlinear_solver) != nullptr)
#endif  // USE_PETSC
    )
    {
        // The Newton-Raphson method needs a Newton-ready ODE.

        using ODENewton = NumLib::ODESystem<ODETag, Tag::Newton>;
        if (auto* ode_newton = dynamic_cast<ODENewton*>(&ode_sys))
        {
            process_data.tdisc_ode_sys = make_tdisc_ode_sys(*ode_newton);
        }
        else
        {
            OGS_FATAL(
                "You are trying to solve a non-Newton-ready ODE with the"
                " Newton-Raphson method. Aborting");
        }
    }
    else
    {
        OGS_FATAL("Encountered unknown nonlinear solver type. Aborting");
    }
}

void setTimeDiscretizedODESystem(ProcessData& process_data,
                                 ProcessData const* const storage_owner)
{
    setTimeDiscretizedODESystem(process_data, process_data.process,
                                storage_owner);
}

void setEquationSystem(ProcessData const& process_data)
{
    auto& conv_crit = *process_data.conv_crit;
    auto& nonlinear_solver = process_data.nonlinear_solver;

    using Tag = NumLib::NonlinearSolverTag;
    switch (process_data.nonlinear_solver_tag)
    {
        case Tag::Picard:
        {
            using EqSys = NumLib::NonlinearSystem<Tag::Picard>;
            auto& eq_sys_ = static_cast<EqSys&>(*process_data.tdisc_ode_sys);
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
            auto& eq_sys_ = static_cast<EqSys&>(*process_data.tdisc_ode_sys);

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
