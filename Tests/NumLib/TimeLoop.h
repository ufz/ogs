#pragma once

#include "TimeDiscretizedODESystem.h"
#include "NonlinearSolver.h"

template<NonlinearSolverTag NLTag>
class TimeLoop
{
public:
    using TDiscODESys = TimeDiscretizedODESystemBase<NLTag>;
    using NLSolver = NonlinearSolver<NLTag>;

    explicit
    TimeLoop(TDiscODESys& ode_sys, NLSolver& nonlinear_solver)
        : _ode_sys(ode_sys)
        , _nonlinear_solver(nonlinear_solver)
    {}

    template<typename Callback>
    void loop(const double t0, const Vector x0,
              const double t_end, const double delta_t,
              Callback& post_timestep);

private:
    TDiscODESys& _ode_sys;
    NLSolver& _nonlinear_solver;
};


template<NonlinearSolverTag NLTag>
template<typename Callback>
void
TimeLoop<NLTag>::
loop(const double t0, const Vector x0, const double t_end, const double delta_t,
     Callback& post_timestep)
{
    Vector x(x0); // solution vector

    auto& time_disc = _ode_sys.getTimeDiscretization();

    time_disc.setInitialState(t0, x0); // push IC

    if (time_disc.needsPreload()) {
        _nonlinear_solver.assemble(_ode_sys, x);
        time_disc.pushState(t0, x0, _ode_sys); // TODO: that might do duplicate work
    }

    for (double t=t0+delta_t; t<t_end+delta_t; t+=delta_t)
    {
        // INFO("time: %e, delta_t: %e", t, delta_t);
        time_disc.setCurrentTime(t, delta_t);

        _nonlinear_solver.solve(_ode_sys, x);

        time_disc.pushState(t, x, _ode_sys);

        post_timestep(t, x);
    }
}
