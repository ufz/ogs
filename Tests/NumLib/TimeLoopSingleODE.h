/**
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "NumLib/DOF/GlobalMatrixProviders.h"
#include "NumLib/ODESolver/TimeDiscretizedODESystem.h"
#include "NumLib/ODESolver/NonlinearSolver.h"
#include "MathLib/LinAlg/LinAlg.h"

namespace NumLib
{
//! \addtogroup ODESolver
//! @{

/*! Integrate a single first-order ODE system over time.
 *
 */
template <NonlinearSolverTag NLTag>
class TimeLoopSingleODE final
{
public:
    using TDiscODESys = TimeDiscretizedODESystemBase<NLTag>;
    using NLSolver = NonlinearSolver<NLTag>;

    /*! Constructs an new instance.
     *
     * \param ode_sys The ODE system to be integrated
     * \param linear_solver the linear solver used to solve the linearized ODE
     *                      system.
     * \param nonlinear_solver The solver to be used to resolve nonlinearities.
     * \param convergence_criterion The convergence criterion used by the
     * nonlinear solver.
     */
    TimeLoopSingleODE(
        TDiscODESys& ode_sys,
        std::unique_ptr<GlobalLinearSolver>&& linear_solver,
        std::unique_ptr<NLSolver>&& nonlinear_solver,
        std::unique_ptr<ConvergenceCriterion>&& convergence_criterion)
        : ode_sys_(ode_sys),
          linear_solver_(std::move(linear_solver)),
          nonlinear_solver_(std::move(nonlinear_solver)),
          convergence_criterion_(std::move(convergence_criterion))
    {
    }

    /*! Integrate the ODE from \c t0 to \c t_end with a timestep size of
     * \c delta_t.
     *
     * The initial condition is \f$ x(\mathtt{t0}) = \mathtt{x0} \f$.
     *
     * After each timestep the callback \c post_timestep will be called,
     *  i.e., it won't be called with the initial condition as parameters.
     *
     * \tparam Callback Any callable object which can be called with the
     *         arguments of type \c double and <tt>Vector const&</tt> which
     *         contain the time and solution at the current timestep.
     *
     * \retval true  if the ODE could be successfully integrated
     * \retval false otherwise
     */
    template <typename Callback>
    NumLib::NonlinearSolverStatus loop(const double t0, GlobalVector const& x0,
                                       const double t_end, const double delta_t,
                                       Callback& post_timestep);

private:
    TDiscODESys& ode_sys_;
    std::unique_ptr<GlobalLinearSolver> linear_solver_;
    std::unique_ptr<NLSolver> nonlinear_solver_;
    std::unique_ptr<ConvergenceCriterion> convergence_criterion_;
};

//! @}

template <NonlinearSolverTag NLTag>
template <typename Callback>
NumLib::NonlinearSolverStatus TimeLoopSingleODE<NLTag>::loop(
    const double t0, GlobalVector const& x0, const double t_end,
    const double delta_t, Callback& post_timestep)
{
    // solution vector
    std::vector<GlobalVector*> xs;
    xs.push_back(&NumLib::GlobalVectorProvider::provider.getVector(x0));
    GlobalVector& x = *xs.back();

    std::vector<GlobalVector*> xs_prev;
    xs_prev.push_back(&NumLib::GlobalVectorProvider::provider.getVector(x0));
    GlobalVector& x_prev = *xs_prev.back();

    auto& time_disc = ode_sys_.getTimeDiscretization();

    time_disc.setInitialState(t0);     // push IC
    MathLib::LinAlg::copy(x, x_prev);  // pushState

    nonlinear_solver_->setEquationSystem(ode_sys_, *convergence_criterion_);

    double t;
    unsigned timestep = 0;
    NumLib::NonlinearSolverStatus nonlinear_solver_status = {false, 0};
    for (t = t0 + delta_t; t < t_end + std::numeric_limits<double>::epsilon();
         t = t0 + (timestep + 1) * delta_t)
    {
        ++timestep;

        // INFO("time: {:e}, delta_t: {:e}", t, delta_t);
        time_disc.nextTimestep(t, delta_t);

        int const process_id = 0;
        nonlinear_solver_status =
            nonlinear_solver_->solve(xs, xs_prev, nullptr, process_id);
        if (!nonlinear_solver_status.error_norms_met)
        {
            break;
        }

        MathLib::LinAlg::copy(x, x_prev);  // pushState

        auto const t_cb =
            t;  // make sure the callback cannot overwrite anything.
        auto const& x_cb = x;  // ditto.
        post_timestep(t_cb, x_cb);
    }

    NumLib::GlobalVectorProvider::provider.releaseVector(x);

    if (!nonlinear_solver_status.error_norms_met)
    {
        ERR("Nonlinear solver failed in timestep #{:d} at t = {:g} s", timestep,
            t);
    }
    return nonlinear_solver_status;
}
}  // namespace NumLib
