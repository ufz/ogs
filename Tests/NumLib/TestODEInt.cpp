#include <gtest/gtest.h>

#include "TimeLoop.h"

#include <logog/include/logog.hpp>


// ODE 1 //////////////////////////////////////////////////////////
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

        Eigen::Matrix2d k;
        k << 0.0, 1.0,
            -1.0, 0.0;

        Jac = (m*dxdot_dx + k).sparseView();
    }

    IndexType getMatrixSize() const override
    {
        return 2;
    }

    bool isLinear() const override
    {
        return true;
    }
};

Eigen::Vector2d ode1_solution(const double t)
{
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
void test_Ode1(TimeDisc& timeDisc)
{
    Ode1 ode;
    TimeDiscretizedODESystem<NLTag, TimeDisc> ode_sys(ode, timeDisc);

    const double tol = 1e-4;
    const unsigned maxiter = 5;
    NonlinearSolver<NLTag> nonlinear_solver(tol, maxiter);

    TimeLoop<NLTag, TimeDisc> loop(ode_sys, nonlinear_solver);

    // initial condition
    const double t0 = 0.0;
    Eigen::Vector2d x0;
    x0 << 1.0, 0.0;

    const double t_end = 0.1;
    const double delta_t = (t_end-t0)/10.0;

    loop.loop(t0, x0, t_end, delta_t, print_result_Ode1);
}
// ODE 1 end //////////////////////////////////////////////////////

// ODE 2 //////////////////////////////////////////////////////////
class Ode2 final
        : public IFirstOrderImplicitOde<NonlinearSolverTag::Newton>
{
public:
    void assemble(const double t, Vector const& x,
                  Matrix& M, Matrix& K, Vector& b) override
    {
        (void) t;
        M.coeffRef(0, 0) = 1.0;
        K.coeffRef(0, 0) = x[0];
        b.coeffRef(0, 0) = 0.0;
    }

    void assembleJacobian(const double t, const Vector &x,
                          const double dxdot_dx, Matrix &Jac)
    {
        (void) t;
        Jac.coeffRef(0, 0) = dxdot_dx + 2.0*x[0];
    }

    IndexType getMatrixSize() const override
    {
        return 1;
    }

    bool isLinear() const override
    {
        return false;
    }
};

double ode2_solution(const double t)
{
    return 1.0/t;
}

void print_result_Ode2(const double t, Vector const& x)
{
    INFO("x[0] = %10g", x[0]);
    auto const v = ode2_solution(t);
    INFO("x[0] = %10g (analytical solution)", v);
}

template<NonlinearSolverTag NLTag, typename TimeDisc>
void test_Ode2(TimeDisc& timeDisc)
{
    Ode2 ode;
    TimeDiscretizedODESystem<NLTag, TimeDisc> ode_sys(ode, timeDisc);

    const double tol = 1e-4;
    const unsigned maxiter = 5;
    NonlinearSolver<NLTag> nonlinear_solver(tol, maxiter);

    TimeLoop<NLTag, TimeDisc> loop(ode_sys, nonlinear_solver);

    // initial condition
    const double t0 = 1.0;
    Eigen::VectorXd x0(1);
    x0[0] = 1.0;

    const double t_end = 2.0;
    const double delta_t = (t_end-t0)/10.0;

    loop.loop(t0, x0, t_end, delta_t, print_result_Ode2);
}
// ODE 2 end //////////////////////////////////////////////////////


TEST(NumLibODEInt, Ode1PicardNewtonBwdEuler)
{
    auto const NLTag = NonlinearSolverTag::Picard;
    using TimeDisc = BackwardEuler;
    TimeDisc timeDisc;
    test_Ode1<NLTag, TimeDisc>(timeDisc);
}


TEST(NumLibODEInt, Ode1PicardNewtonFwdEuler)
{
    auto const NLTag = NonlinearSolverTag::Newton;
    using TimeDisc = ForwardEuler;
    TimeDisc timeDisc;
    test_Ode1<NLTag, TimeDisc>(timeDisc);
}


TEST(NumLibODEInt, Ode2PicardNewtonBwdEuler)
{
    auto const NLTag = NonlinearSolverTag::Newton;
    using TimeDisc = BackwardEuler;
    TimeDisc timeDisc;
    test_Ode2<NLTag, TimeDisc>(timeDisc);
}


TEST(NumLibODEInt, Ode2PicardNewtonCrankNicolson)
{
    auto const NLTag = NonlinearSolverTag::Newton;
    using TimeDisc = CrankNicolson;
    TimeDisc timeDisc;
    test_Ode2<NLTag, TimeDisc>(timeDisc);
}


