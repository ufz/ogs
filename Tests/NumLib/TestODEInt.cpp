#include <gtest/gtest.h>
#include <logog/include/logog.hpp>

#include <fstream>
#include <memory>
#include <typeinfo>

#include "TimeLoop.h"
#include "Odes.h"
#include "BaseLib/BuildInfo.h"


template<NonlinearSolverTag NLTag>
class TestOutput
{
public:
    template<typename Ode>
    void run_test(Ode& ode, TimeDiscretization& timeDisc)
    {
        run_test<Ode>(ode, timeDisc, 10); // by default make 10 timesteps
    }

    template<typename Ode>
    void run_test(Ode& ode, TimeDiscretization& timeDisc, const unsigned num_timesteps)
    {
        auto mat_trans = createMatrixTranslator<ODESystemTag::FirstOrderImplicitQuasilinear>(timeDisc);
        TimeDiscretizedODESystem<ODESystemTag::FirstOrderImplicitQuasilinear, NLTag> ode_sys(ode, timeDisc, *mat_trans);
        TimeLoop<NLTag> loop(ode_sys, _nonlinear_solver);

        const double t0      = OdeTraits<Ode>::t0;
        const double t_end   = OdeTraits<Ode>::t_end;
        const double delta_t = (t_end-t0) / num_timesteps;

        init_file(ode, timeDisc, delta_t);

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
    template<typename Ode, typename TimeDisc>
    void init_file(Ode const& ode, TimeDisc const& timeDisc, const double delta_t)
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
        for (IndexType i=0; i<x_num.size(); ++i) *_file << '\t' << x_num[i];
        for (IndexType i=0; i<x_ana.size(); ++i) *_file << '\t' << x_ana[i];
        *_file << "\n";
    }

    template<typename Ode>
    void loopCallback(const double t, Vector const& x)
    {
        write(t, x, OdeTraits<Ode>::solution(t));
    }

    std::unique_ptr<std::ofstream> _file;

    const double _tol = 1e-8;
    const unsigned _maxiter = 10;
    NonlinearSolver<NLTag> _nonlinear_solver = NonlinearSolver<NLTag>(_tol, _maxiter);
};


template<typename TimeDisc, typename ODE, NonlinearSolverTag NLTag>
typename std::enable_if<std::is_default_constructible<TimeDisc>::value>::type
run_test_case(const unsigned num_timesteps)
{
    ODE ode;
    TimeDisc timeDisc;

    TestOutput<NLTag> test;
    test.run_test(ode, timeDisc, num_timesteps);
}

template<typename TimeDisc, typename ODE, NonlinearSolverTag NLTag>
typename std::enable_if<std::is_same<TimeDisc, CrankNicolson>::value>::type
run_test_case(const unsigned num_timesteps)
{
    ODE ode;
    TimeDisc timeDisc(0.5);

    TestOutput<NLTag> test;
    test.run_test(ode, timeDisc, num_timesteps);
}

template<typename TimeDisc, typename ODE, NonlinearSolverTag NLTag>
typename std::enable_if<std::is_same<TimeDisc, BackwardDifferentiationFormula>::value>::type
run_test_case(const unsigned num_timesteps)
{
    ODE ode;
    TimeDisc timeDisc(3);

    TestOutput<NLTag> test;
    test.run_test(ode, timeDisc, num_timesteps);
}



template<typename ODE_, typename TimeDisc_, NonlinearSolverTag NLTag_>
struct TestCase
{
    using ODE = ODE_;
    using TimeDisc = TimeDisc_;
    static constexpr NonlinearSolverTag NLTag = NLTag_;
};


typedef ::testing::Types<
    TestCase<Ode1, BackwardEuler,                  NonlinearSolverTag::Newton>,
    TestCase<Ode1, ForwardEuler,                   NonlinearSolverTag::Newton>,
    TestCase<Ode1, CrankNicolson,                  NonlinearSolverTag::Newton>,
    TestCase<Ode1, BackwardDifferentiationFormula, NonlinearSolverTag::Newton>,

    TestCase<Ode1, BackwardEuler,                  NonlinearSolverTag::Picard>,
    TestCase<Ode1, ForwardEuler,                   NonlinearSolverTag::Picard>,
    TestCase<Ode1, CrankNicolson,                  NonlinearSolverTag::Picard>,
    TestCase<Ode1, BackwardDifferentiationFormula, NonlinearSolverTag::Picard>,

    TestCase<Ode2, BackwardEuler,                  NonlinearSolverTag::Newton>,
    TestCase<Ode2, ForwardEuler,                   NonlinearSolverTag::Newton>,
    TestCase<Ode2, CrankNicolson,                  NonlinearSolverTag::Newton>,
    TestCase<Ode2, BackwardDifferentiationFormula, NonlinearSolverTag::Newton>,

    TestCase<Ode2, BackwardEuler,                  NonlinearSolverTag::Picard>,
    TestCase<Ode2, ForwardEuler,                   NonlinearSolverTag::Picard>,
    TestCase<Ode2, CrankNicolson,                  NonlinearSolverTag::Picard>,
    TestCase<Ode2, BackwardDifferentiationFormula, NonlinearSolverTag::Picard>
> TestCases;



template<class TestParams>
class NumLibODEIntTyped : public ::testing::Test
{
public:
    using ODE      = typename TestParams::ODE;
    using TimeDisc = typename TestParams::TimeDisc;
    static constexpr NonlinearSolverTag NLTag = TestParams::NLTag;

    static void test()
    {
        for (auto num_timesteps : { 10, 100, 1000 }) {
            run_test_case<TimeDisc, ODE, NLTag>(num_timesteps);
        }
    }
};

TYPED_TEST_CASE(NumLibODEIntTyped, TestCases);

TYPED_TEST(NumLibODEIntTyped, T1)
{
    TestFixture::test();
}
