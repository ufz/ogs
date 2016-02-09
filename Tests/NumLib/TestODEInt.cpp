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
    static void setIC(Vector& x0)
    {
        x0.resize(2);
        x0 << 1.0, 0.0;
    }

    static Vector solution(const double t)
    {
        Eigen::Vector2d v;
        v << cos(t), sin(t);
        return v;
    }

    static constexpr double t0    = 0.0;
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

template<>
class OdeTraits<Ode2>
{
public:
    static void setIC(Vector& x0)
    {
        x0.resize(1);
        x0[0] = 1.0;
    }

    static Vector solution(const double t)
    {
        Vector v(1);
        v[0] = 1.0 / t;
        return v;
    }

    static constexpr double t0    = 1.0;
    static constexpr double t_end = 2.0;
};
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
        auto const t_end = OdeTraits<Ode>::t_end;
        auto const t0    = OdeTraits<Ode>::t0;
        run_test<Ode>(ode, timeDisc, (t_end-t0)/10.0); // by default make 10 timesteps
    }

    template<typename Ode>
    void run_test(Ode& ode, ITimeDiscretization& timeDisc, const double delta_t)
    {
        TimeDiscretizedODESystem<NLTag> ode_sys(ode, timeDisc);
        TimeLoop<NLTag> loop(ode_sys, _nonlinear_solver);

        const double t0    = OdeTraits<Ode>::t0;
        const double t_end = OdeTraits<Ode>::t_end;

        // initial condition
        Vector x0;
        OdeTraits<Ode>::setIC(x0);

        write(t0, x0, x0);

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

    const double _tol = 1e-8;
    const unsigned _maxiter = 10;
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

    Ode2 ode;
    TimeDisc timeDisc;

    TestOutput<NLTag> test("Ode2_BwdEuler");
    test.run_test(ode, timeDisc);
}


TEST(NumLibODEInt, Ode2_CrankNicolson)
{
    auto const NLTag = NonlinearSolverTag::Newton;
    using TimeDisc = CrankNicolson;

    Ode2 ode;
    TimeDisc timeDisc(0.5);

    TestOutput<NLTag> test("Ode2_CrankNicolson");
    test.run_test(ode, timeDisc);
}


TEST(NumLibODEInt, Ode2_BDF)
{
    auto const NLTag = NonlinearSolverTag::Newton;
    using TimeDisc = BackwardDifferentiationFormula;

    Ode2 ode;
    TimeDisc timeDisc(3);

    TestOutput<NLTag> test("Ode2_BDF");
    test.run_test(ode, timeDisc);
}


