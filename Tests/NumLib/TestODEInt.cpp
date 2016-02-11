#include <gtest/gtest.h>
#include <logog/include/logog.hpp>

#include <fstream>
#include <memory>
#include <typeinfo>

#include "TimeLoop.h"
#include "Odes.h"
#include "BaseLib/BuildInfo.h"


template<typename Matrix, typename Vector, NonlinearSolverTag NLTag>
class TestOutput
{
public:
    using TimeDisc = TimeDiscretization<Vector>;

    template<typename ODE>
    void run_test(ODE& ode, TimeDisc& timeDisc)
    {
        run_test<ODE>(ode, timeDisc, 10); // by default make 10 timesteps
    }

    template<template<typename /*Matrix*/, typename /*Vector*/> typename ODE>
    void run_test(ODE<Matrix, Vector>& ode, TimeDisc& timeDisc, const unsigned num_timesteps)
    {
        using ODE_ = ODE<Matrix, Vector>;
        auto mat_trans = createMatrixTranslator<Matrix, Vector, ODE_::ODETag>(timeDisc);
        TimeDiscretizedODESystem<Matrix, Vector, ODE_::ODETag, NLTag>
                ode_sys(ode, timeDisc, *mat_trans);
        TimeLoop<Matrix, Vector, NLTag> loop(ode_sys, _nonlinear_solver);

        const double t0      = OdeTraits<Vector, ODE>::t0;
        const double t_end   = OdeTraits<Vector, ODE>::t_end;
        const double delta_t = (t_end-t0) / num_timesteps;

        init_file(ode, timeDisc, delta_t);

        // initial condition
        Vector x0(ode.getMatrixSize());
        OdeTraits<Vector, ODE>::setIC(x0);

        write(t0, x0, x0);

        auto cb = [this](const double t, Vector const& x) {
            loopCallback<ODE>(t, x);
        };
        loop.loop(t0, x0, t_end, delta_t, cb);
    }

private:
    template<typename ODE, typename TimeDisc_>
    void init_file(ODE const& ode, TimeDisc_ const& timeDisc, const double delta_t)
    {
        std::string path(BaseLib::BuildInfo::tests_tmp_path + "ODEInt_");
        path += typeid(ode).name();
        path += "_";
        path += typeid(timeDisc).name();
        path += "_";

        switch (NLTag) {
        case NonlinearSolverTag::Picard: path += "Picard"; break;
        case NonlinearSolverTag::Newton: path += "Newton"; break;
        }

        path += "_" + std::to_string(delta_t);
        path += ".csv";

        _file.reset(new std::ofstream(path));
        _file->precision(15);
    }

    void write(double const t, Vector const& x_num, Vector const& x_ana)
    {
        *_file << t;
        for (decltype(x_num.size()) i=0; i<x_num.size(); ++i) *_file << '\t' << x_num[i];
        for (decltype(x_num.size()) i=0; i<x_ana.size(); ++i) *_file << '\t' << x_ana[i];
        *_file << "\n";
    }

    template<template<typename /*Matrix*/, typename /*Vector*/> typename Ode>
    void loopCallback(const double t, Vector const& x)
    {
        write(t, x, OdeTraits<Vector, Ode>::solution(t));
    }

    std::unique_ptr<std::ofstream> _file;

    const double _tol = 1e-8;
    const unsigned _maxiter = 10;

    using NLSolver = NonlinearSolver<Matrix, Vector, NLTag>;
    NLSolver _nonlinear_solver = NLSolver(_tol, _maxiter);
};


template<typename Matrix, typename Vector, typename TimeDisc, typename ODE, NonlinearSolverTag NLTag>
typename std::enable_if<std::is_default_constructible<TimeDisc>::value>::type
run_test_case(const unsigned num_timesteps)
{
    ODE ode;
    TimeDisc timeDisc;

    TestOutput<Matrix, Vector, NLTag> test;
    test.run_test(ode, timeDisc, num_timesteps);
}

template<typename Matrix, typename Vector, typename TimeDisc, typename ODE, NonlinearSolverTag NLTag>
typename std::enable_if<std::is_same<TimeDisc, CrankNicolson<Vector> >::value>::type
run_test_case(const unsigned num_timesteps)
{
    ODE ode;
    TimeDisc timeDisc(0.5);

    TestOutput<Matrix, Vector, NLTag> test;
    test.run_test(ode, timeDisc, num_timesteps);
}

template<typename Matrix, typename Vector, typename TimeDisc, typename ODE, NonlinearSolverTag NLTag>
typename std::enable_if<
    std::is_same<TimeDisc, BackwardDifferentiationFormula<Vector> >::value>::type
run_test_case(const unsigned num_timesteps)
{
    ODE ode;
    TimeDisc timeDisc(3);

    TestOutput<Matrix, Vector, NLTag> test;
    test.run_test(ode, timeDisc, num_timesteps);
}



template<typename Matrix_, typename Vector_,
         template<typename /*Matrix*/, typename /*Vector*/> typename ODE_,
         template<typename /*Vector*/> typename TimeDisc_,
         NonlinearSolverTag NLTag_>
struct TestCase
{
    using Matrix = Matrix_;
    using Vector = Vector_;
    using ODE = ODE_<Matrix_, Vector_>;
    using TimeDisc = TimeDisc_<Vector_>;
    static constexpr NonlinearSolverTag NLTag = NLTag_;
};


typedef ::testing::Types<
    TestCase<ODEMatrix, ODEVector, Ode1, BackwardEuler,                  NonlinearSolverTag::Newton>,
    TestCase<ODEMatrix, ODEVector, Ode1, ForwardEuler,                   NonlinearSolverTag::Newton>,
    TestCase<ODEMatrix, ODEVector, Ode1, CrankNicolson,                  NonlinearSolverTag::Newton>,
    TestCase<ODEMatrix, ODEVector, Ode1, BackwardDifferentiationFormula, NonlinearSolverTag::Newton>,

    TestCase<ODEMatrix, ODEVector, Ode1, BackwardEuler,                  NonlinearSolverTag::Picard>,
    TestCase<ODEMatrix, ODEVector, Ode1, ForwardEuler,                   NonlinearSolverTag::Picard>,
    TestCase<ODEMatrix, ODEVector, Ode1, CrankNicolson,                  NonlinearSolverTag::Picard>,
    TestCase<ODEMatrix, ODEVector, Ode1, BackwardDifferentiationFormula, NonlinearSolverTag::Picard>,

    TestCase<ODEMatrix, ODEVector, Ode2, BackwardEuler,                  NonlinearSolverTag::Newton>,
    TestCase<ODEMatrix, ODEVector, Ode2, ForwardEuler,                   NonlinearSolverTag::Newton>,
    TestCase<ODEMatrix, ODEVector, Ode2, CrankNicolson,                  NonlinearSolverTag::Newton>,
    TestCase<ODEMatrix, ODEVector, Ode2, BackwardDifferentiationFormula, NonlinearSolverTag::Newton>,

    TestCase<ODEMatrix, ODEVector, Ode2, BackwardEuler,                  NonlinearSolverTag::Picard>,
    TestCase<ODEMatrix, ODEVector, Ode2, ForwardEuler,                   NonlinearSolverTag::Picard>,
    TestCase<ODEMatrix, ODEVector, Ode2, CrankNicolson,                  NonlinearSolverTag::Picard>,
    TestCase<ODEMatrix, ODEVector, Ode2, BackwardDifferentiationFormula, NonlinearSolverTag::Picard>
> TestCases;



template<class TestParams>
class NumLibODEIntTyped : public ::testing::Test
{
public:
    using Matrix   = typename TestParams::Matrix;
    using Vector   = typename TestParams::Vector;
    using ODE      = typename TestParams::ODE;
    using TimeDisc = typename TestParams::TimeDisc;
    static constexpr NonlinearSolverTag NLTag = TestParams::NLTag;

    static void test()
    {
        for (auto num_timesteps : { 10, 100, 1000 }) {
            run_test_case<Matrix, Vector, TimeDisc, ODE, NLTag>(num_timesteps);
        }
    }
};

TYPED_TEST_CASE(NumLibODEIntTyped, TestCases);

TYPED_TEST(NumLibODEIntTyped, T1)
{
    TestFixture::test();
}
