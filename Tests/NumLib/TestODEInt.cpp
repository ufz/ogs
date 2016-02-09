#include <gtest/gtest.h>
#include <logog/include/logog.hpp>

// #include <iomanip>
#include <fstream>

#include "TimeLoop.h"
#include "BaseLib/BuildInfo.h"

template<typename Ode>
class OdeTraits;

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

template<>
class OdeTraits<Ode1>
{
public:
    static void setIC(double& t0, Vector& x0)
    {
        t0 = 0.0;
        x0.resize(2);
        x0 << 1.0, 0.0;
    }

    static Vector solution(const double t)
    {
        Eigen::Vector2d v;
        v << cos(t), sin(t);
        return v;
    }

    static constexpr double t_end = 1.0;
};
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

template<NonlinearSolverTag NLTag>
void test_Ode2(ITimeDiscretization& timeDisc)
{
    Ode2 ode;
    TimeDiscretizedODESystem<NLTag> ode_sys(ode, timeDisc);

    const double tol = 1e-4;
    const unsigned maxiter = 5;
    NonlinearSolver<NLTag> nonlinear_solver(tol, maxiter);

    TimeLoop<NLTag> loop(ode_sys, nonlinear_solver);

    // initial condition
    const double t0 = 1.0;
    Eigen::VectorXd x0(1);
    x0[0] = 1.0;

    const double t_end = 2.0;
    const double delta_t = (t_end-t0)/10.0;

    loop.loop(t0, x0, t_end, delta_t, print_result_Ode2);
}
// ODE 2 end //////////////////////////////////////////////////////



template<NonlinearSolverTag NLTag>
class TestOutput
{
public:
    TestOutput(std::string const& filename)
        : _file(BaseLib::BuildInfo::tests_tmp_path + "ODEInt_" + filename + ".csv")
    {
        _file.precision(15);
    }

    template<typename Ode>
    void run_test(Ode& ode, ITimeDiscretization& timeDisc)
    {
        TimeDiscretizedODESystem<NLTag> ode_sys(ode, timeDisc);
        TimeLoop<NLTag> loop(ode_sys, _nonlinear_solver);

        // initial condition
        double t0;
        Vector x0;
        OdeTraits<Ode>::setIC(t0, x0);

        write(t0, x0, x0);

        const double t_end = OdeTraits<Ode>::t_end;
        const double delta_t = (t_end-t0)/10.0;

        auto cb = [this](const double t, Vector const& x) {
            loopCallback<Ode>(t, x);
        };
        loop.loop(t0, x0, t_end, delta_t, cb);
    }

private:
    void write(double const t, Vector const& x_num, Vector const& x_ana)
    {
        _file << t;
        for (IndexType i=0; i<x_ana.size(); ++i) _file << '\t' << x_ana[i];
        for (IndexType i=0; i<x_num.size(); ++i) _file << '\t' << x_num[i];
        _file << "\n";
    }

    template<typename Ode>
    void loopCallback(const double t, Vector const& x)
    {
        write(t, x, OdeTraits<Ode>::solution(t));
    }


    std::ofstream _file;

    const double _tol = 1e-4;
    const unsigned _maxiter = 5;
    NonlinearSolver<NLTag> _nonlinear_solver = NonlinearSolver<NLTag>(_tol, _maxiter);
};


TEST(NumLibODEInt, Ode1_BwdEuler)
{
    auto const NLTag = NonlinearSolverTag::Picard;
    using TimeDisc = BackwardEuler;

    Ode1 ode;
    TimeDisc timeDisc;

    TestOutput<NLTag> test("Ode1_BwdEuler");
    test.run_test(ode, timeDisc);
}


TEST(NumLibODEInt, Ode1_FwdEuler)
{
    auto const NLTag = NonlinearSolverTag::Newton;
    using TimeDisc = ForwardEuler;

    Ode1 ode;
    TimeDisc timeDisc;

    TestOutput<NLTag> test("Ode1_FwdEuler");
    test.run_test(ode, timeDisc);
}


TEST(NumLibODEInt, Ode2_BwdEuler)
{
    auto const NLTag = NonlinearSolverTag::Newton;
    using TimeDisc = BackwardEuler;
    TimeDisc timeDisc;
    test_Ode2<NLTag>(timeDisc);
}


TEST(NumLibODEInt, Ode2_CrankNicolson)
{
    auto const NLTag = NonlinearSolverTag::Newton;
    using TimeDisc = CrankNicolson;
    TimeDisc timeDisc(0.5);
    test_Ode2<NLTag>(timeDisc);
}


TEST(NumLibODEInt, Ode2_BDF)
{
    auto const NLTag = NonlinearSolverTag::Newton;
    using TimeDisc = BackwardDifferentiationFormula;
    TimeDisc timeDisc(3);
    test_Ode2<NLTag>(timeDisc);
}


