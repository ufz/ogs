#include <gtest/gtest.h>

#include "NonlinSolver.h"
#include <logog/include/logog.hpp>


class Ode1 final
        : public IFirstOrderImplicitOde<NonlinearSolverTag::Newton>
{
public:
    void assemble(const double t, Vector const& x,
                  Matrix& M, Matrix& K, Vector& b) override
    {
        (void) t; (void) x;

        Eigen::Matrix2d m;
        m << 1.0, 0.0,
             0.0, 1.0;
        M = m.sparseView();

        Eigen::Matrix2d k;
        k << 0.0, 1.0,
            -1.0, 0.0;
        K = k.sparseView();

        Eigen::Vector2d b_;
        b_ << 0.0, 0.0;
        b = b_;
    }

    void assembleJacobian(const double t, const Vector &x,
                          const double dxdot_dx, Matrix &Jac)
    {
        (void) t; (void) x;

        Eigen::Matrix2d m;
        m << 1.0, 0.0,
             0.0, 1.0;

        // TODO: \dot x contribution
        Eigen::Matrix2d k;
        k << 0.0, 1.0,
            -1.0, 0.0;

        Jac = (m*dxdot_dx + k).sparseView();
    }

    bool isLinear() const override
    {
        return true;
    }
};

Eigen::Vector2d ode1_solution(const double t) {
    Eigen::Vector2d v;
    v << cos(t), sin(t);
    return v;
}

void print_result_Ode1(const double t, Vector const& x)
{
    INFO("x[0] = %10g, x[1] = %10g", x[0], x[1]);
    auto const v = ode1_solution(t);
    INFO("x[0] = %10g, x[1] = %10g (analytical solution)", v[0], v[1]);
}

template<NonlinearSolverTag NLTag, typename TimeDisc>
void test_Ode1()
{
    Ode1 ode;
    TimeDiscretizedODESystem<NLTag, TimeDisc> ode_sys(ode);

    const double tol = 1e-4;
    const unsigned maxiter = 5;
    NonlinearSolver<NLTag> nonlinear_solver(tol, maxiter);

    TimeLoop<NLTag, TimeDisc> loop(ode_sys, nonlinear_solver);

    const double t_end = 0.1;
    const double delta_t = t_end/10.0;

    // initial condition
    const double t0 = 0.0;
    Eigen::Vector2d x0;
    x0 << 1.0, 0.0;

    loop.loop(t0, x0, t_end, delta_t, print_result_Ode1);
}


TEST(NumLibODEInt, PicardNewtonBwdEuler)
{
    auto const NLTag = NonlinearSolverTag::Newton;
    using TimeDisc = BackwardEuler;
    test_Ode1<NLTag, TimeDisc>();
}


TEST(NumLibODEInt, PicardNewtonFwdEuler)
{
    auto const NLTag = NonlinearSolverTag::Newton;
    using TimeDisc = ForwardEuler;
    test_Ode1<NLTag, TimeDisc>();
}


