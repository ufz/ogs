#pragma once

#include "TimeDiscretizedODESystem.h"
#include "NonlinearSolver.h"


namespace NumLib
{

template<typename Matrix, typename Vector, NonlinearSolverTag NLTag>
class TimeLoop
{
public:
    using TDiscODESys = TimeDiscretizedODESystemBase<Matrix, Vector, NLTag>;
    using NLSolver = NonlinearSolver<Matrix, Vector, NLTag>;

    explicit
    TimeLoop(TDiscODESys& ode_sys, NLSolver& nonlinear_solver)
        : _ode_sys(ode_sys)
        , _nonlinear_solver(nonlinear_solver)
    {}

    template<typename Callback>
    bool loop(const double t0, const Vector x0,
              const double t_end, const double delta_t,
              Callback& post_timestep);

private:
    TDiscODESys& _ode_sys;
    NLSolver& _nonlinear_solver;
};


template<typename Matrix, typename Vector, NonlinearSolverTag NLTag>
template<typename Callback>
bool
TimeLoop<Matrix, Vector, NLTag>::
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

    double t;
    unsigned timestep = 0;
    bool nl_slv_succeeded = true;
    for (t=t0+delta_t; t<t_end+delta_t; t+=delta_t, ++timestep)
    {
        // INFO("time: %e, delta_t: %e", t, delta_t);
        time_disc.setCurrentTime(t, delta_t);

        nl_slv_succeeded = _nonlinear_solver.solve(_ode_sys, x);
        if (!nl_slv_succeeded) break;

        time_disc.pushState(t, x, _ode_sys);

        post_timestep(t, x);
    }

    if (!nl_slv_succeeded) {
        ERR("Nonlinear solver failed in timestep #%u at t = %g s", timestep, t);
    }
    return nl_slv_succeeded;
}

}
