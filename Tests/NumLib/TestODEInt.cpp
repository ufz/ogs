#include <gtest/gtest.h>
#include <logog/include/logog.hpp>

#include <fstream>
#include <memory>
#include <typeinfo>

#include "ProcessLib/NumericsConfig.h"
#include "NumLib/ODESolver/TimeLoop.h"
#include "ODEs.h"
#include "BaseLib/BuildInfo.h"
#include "NumLib/ODESolver/ODETypes.h"


using EDMatrix = Eigen::MatrixXd;
using EVector = Eigen::VectorXd;

using GMatrix = GlobalSetupType::MatrixType;
using GVector = GlobalSetupType::VectorType;


template<typename Matrix, typename Vector, NumLib::NonlinearSolverTag NLTag>
class TestOutput
{
public:
    using TimeDisc = NumLib::TimeDiscretization<Vector>;

    TestOutput(const char* name)
        : _file_name_part(name)
    {}

    template<typename ODE>
    void run_test(ODE& ode, TimeDisc& timeDisc)
    {
        run_test<ODE>(ode, timeDisc, 10); // by default make 10 timesteps
    }

    template<template<typename /*Matrix*/, typename /*Vector*/> class ODE>
    void run_test(ODE<Matrix, Vector>& ode, TimeDisc& timeDisc, const unsigned num_timesteps)
    {
        using ODE_ = ODE<Matrix, Vector>;
        using ODET = ODETraits<Matrix, Vector, ODE>;

        NumLib::TimeDiscretizedODESystem<Matrix, Vector, ODE_::ODETag, NLTag>
                ode_sys(ode, timeDisc);
        NumLib::TimeLoop<Matrix, Vector, NLTag> loop(ode_sys, _nonlinear_solver);

        const double t0      = ODET::t0;
        const double t_end   = ODET::t_end;
        const double delta_t = (t_end-t0) / num_timesteps;

        INFO("Running test %s with %u timesteps of size %g s.",
             _file_name_part.c_str(), num_timesteps, delta_t);
        init_file(delta_t);

        // initial condition
        Vector x0(ode.getNumEquations());
        ODET::setIC(x0);

        write(t0, x0, x0);

        auto cb = [this](const double t, Vector const& x) {
            loopCallback<ODE>(t, x);
        };

        EXPECT_TRUE(loop.loop(t0, x0, t_end, delta_t, cb));
    }

private:
    void init_file(const double delta_t)
    {
        std::string path(BaseLib::BuildInfo::tests_tmp_path + "ODEInt_");
        path += _file_name_part;
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

    template<template<typename /*Matrix*/, typename /*Vector*/> class Ode>
    void loopCallback(const double t, Vector const& x)
    {
        write(t, x, ODETraits<Matrix, Vector, Ode>::solution(t));
    }

    const std::string _file_name_part;
    std::unique_ptr<std::ofstream> _file;

    const double _tol = 1e-9;
    const unsigned _maxiter = 20;

    NumLib::NonlinearSolver<Matrix, Vector, NLTag> _nonlinear_solver{_tol, _maxiter};
};


// TODO remove debugging macro
#if 1

template<typename Matrix, typename Vector, typename TimeDisc, typename ODE,
         NumLib::NonlinearSolverTag NLTag>
typename std::enable_if<std::is_default_constructible<TimeDisc>::value>::type
run_test_case(const unsigned num_timesteps, const char* name)
{
    ODE ode;
    TimeDisc timeDisc;

    TestOutput<Matrix, Vector, NLTag> test(name);
    test.run_test(ode, timeDisc, num_timesteps);
}

template<typename Matrix, typename Vector, typename TimeDisc, typename ODE,
         NumLib::NonlinearSolverTag NLTag>
typename std::enable_if<std::is_same<TimeDisc, NumLib::CrankNicolson<Vector> >::value>::type
run_test_case(const unsigned num_timesteps, const char* name)
{
    ODE ode;
    TimeDisc timeDisc(0.5);

    TestOutput<Matrix, Vector, NLTag> test(name);
    test.run_test(ode, timeDisc, num_timesteps);
}

template<typename Matrix, typename Vector, typename TimeDisc, typename ODE,
         NumLib::NonlinearSolverTag NLTag>
typename std::enable_if<
    std::is_same<TimeDisc, NumLib::BackwardDifferentiationFormula<Vector> >::value>::type
run_test_case(const unsigned num_timesteps, const char* name)
{
    ODE ode;
    TimeDisc timeDisc(3);

    TestOutput<Matrix, Vector, NLTag> test(name);
    test.run_test(ode, timeDisc, num_timesteps);
}

#else

template<typename Matrix, typename Vector, typename TimeDisc, typename ODE,
         NumLib::NonlinearSolverTag NLTag>
void
run_test_case(const unsigned num_timesteps, const char* name)
{
    ODE ode;
    // TimeDisc timeDisc(3);

    TestOutput<Matrix, Vector, NLTag> test(name);
    // test.run_test(ode, timeDisc, num_timesteps);
}
#endif


// This class is only here s.t. I don't have to put the members into
// the definition of the macro TCLITEM below.
template<typename Matrix_, typename Vector_,
         template<typename /*Matrix*/, typename /*Vector*/> class ODE_,
         template<typename /*Vector*/> class TimeDisc_,
         NumLib::NonlinearSolverTag NLTag_>
struct TestCaseBase
{
    using Matrix = Matrix_;
    using Vector = Vector_;
    using ODE = ODE_<Matrix_, Vector_>;
    using TimeDisc = TimeDisc_<Vector_>;
    static const NumLib::NonlinearSolverTag NLTag = NLTag_;
};


template<typename Matrix_, typename Vector_,
         template<typename /*Matrix*/, typename /*Vector*/> class ODE_,
         template<typename /*Vector*/> class TimeDisc_,
         NumLib::NonlinearSolverTag NLTag_>
struct TestCase;


// /////////////////////////////////////
//
//  Put new test cases to that list
//
// /////////////////////////////////////
#define TESTCASESLIST \
    /* Eigen dense matrix */ \
    TCLITEM(EDMatrix, EVector, ODE1, BackwardEuler,                  Newton) TCLSEP \
    TCLITEM(EDMatrix, EVector, ODE1, ForwardEuler,                   Newton) TCLSEP \
    TCLITEM(EDMatrix, EVector, ODE1, CrankNicolson,                  Newton) TCLSEP \
    TCLITEM(EDMatrix, EVector, ODE1, BackwardDifferentiationFormula, Newton) TCLSEP \
    \
    TCLITEM(EDMatrix, EVector, ODE1, BackwardEuler,                  Picard) TCLSEP \
    TCLITEM(EDMatrix, EVector, ODE1, ForwardEuler,                   Picard) TCLSEP \
    TCLITEM(EDMatrix, EVector, ODE1, CrankNicolson,                  Picard) TCLSEP \
    TCLITEM(EDMatrix, EVector, ODE1, BackwardDifferentiationFormula, Picard) TCLSEP \
    \
    TCLITEM(EDMatrix, EVector, ODE2, BackwardEuler,                  Newton) TCLSEP \
    TCLITEM(EDMatrix, EVector, ODE2, ForwardEuler,                   Newton) TCLSEP \
    TCLITEM(EDMatrix, EVector, ODE2, CrankNicolson,                  Newton) TCLSEP \
    TCLITEM(EDMatrix, EVector, ODE2, BackwardDifferentiationFormula, Newton) TCLSEP \
    \
    TCLITEM(EDMatrix, EVector, ODE2, BackwardEuler,                  Picard) TCLSEP \
    TCLITEM(EDMatrix, EVector, ODE2, ForwardEuler,                   Picard) TCLSEP \
    TCLITEM(EDMatrix, EVector, ODE2, CrankNicolson,                  Picard) TCLSEP \
    TCLITEM(EDMatrix, EVector, ODE2, BackwardDifferentiationFormula, Picard) TCLSEP \
    \
    /* Global sparse matrix */ \
    TCLITEM(GMatrix,  GVector, ODE1, BackwardEuler,                  Newton) TCLSEP \
    TCLITEM(GMatrix,  GVector, ODE1, ForwardEuler,                   Newton) TCLSEP \
    TCLITEM(GMatrix,  GVector, ODE1, CrankNicolson,                  Newton) TCLSEP \
    TCLITEM(GMatrix,  GVector, ODE1, BackwardDifferentiationFormula, Newton) TCLSEP \
    \
    TCLITEM(GMatrix,  GVector, ODE1, BackwardEuler,                  Picard) TCLSEP \
    TCLITEM(GMatrix,  GVector, ODE1, ForwardEuler,                   Picard) TCLSEP \
    TCLITEM(GMatrix,  GVector, ODE1, CrankNicolson,                  Picard) TCLSEP \
    TCLITEM(GMatrix,  GVector, ODE1, BackwardDifferentiationFormula, Picard) TCLSEP \
    \
    TCLITEM(GMatrix,  GVector, ODE2, BackwardEuler,                  Newton) TCLSEP \
    TCLITEM(GMatrix,  GVector, ODE2, ForwardEuler,                   Newton) TCLSEP \
    TCLITEM(GMatrix,  GVector, ODE2, CrankNicolson,                  Newton) TCLSEP \
    TCLITEM(GMatrix,  GVector, ODE2, BackwardDifferentiationFormula, Newton) TCLSEP \
    \
    TCLITEM(GMatrix,  GVector, ODE2, BackwardEuler,                  Picard) TCLSEP \
    TCLITEM(GMatrix,  GVector, ODE2, ForwardEuler,                   Picard) TCLSEP \
    TCLITEM(GMatrix,  GVector, ODE2, CrankNicolson,                  Picard) TCLSEP \
    TCLITEM(GMatrix,  GVector, ODE2, BackwardDifferentiationFormula, Picard)

// ODE3 behaves too badly. Even with very tiny timesteps it cannot be solved.
// Probably because then the singular M matrix has too much weight.
//  TCLITEM(EDMatrix, EVector, ODE3, BackwardEuler,                  Newton) TCLSEP
//  /* Not possible because of singular matrix */
//  /* TCLITEM(EDMatrix, EVector, ODE3, ForwardEuler,                   Newton) TCLSEP */
//  TCLITEM(EDMatrix, EVector, ODE3, CrankNicolson,                  Newton) TCLSEP
//  TCLITEM(EDMatrix, EVector, ODE3, BackwardDifferentiationFormula, Newton) TCLSEP
//
//  TCLITEM(EDMatrix, EVector, ODE3, BackwardEuler,                  Picard) TCLSEP
//  /* Not possible because of singular matrix */
//  /* TCLITEM(EDMatrix, EVector, ODE3, ForwardEuler,                   Picard) TCLSEP */
//  TCLITEM(EDMatrix, EVector, ODE3, CrankNicolson,                  Picard) TCLSEP
//  TCLITEM(EDMatrix, EVector, ODE3, BackwardDifferentiationFormula, Picard)


#define TCLITEM(MAT, VEC, ODE, TIMEDISC, NLTAG) \
    template<> \
    struct TestCase<MAT, VEC, ODE, NumLib::TIMEDISC, NumLib::NonlinearSolverTag::NLTAG> \
        : TestCaseBase<MAT, VEC, ODE, NumLib::TIMEDISC, NumLib::NonlinearSolverTag::NLTAG> \
    { \
        static const char name[]; \
    }; \
    const char TestCase<MAT, VEC, ODE, NumLib::TIMEDISC, \
                        NumLib::NonlinearSolverTag::NLTAG>::name[] \
        = #MAT "_" #VEC "_" #ODE "_" #TIMEDISC "_" #NLTAG;
#define TCLSEP

TESTCASESLIST

#undef TCLITEM
#undef TCLSEP

#define TCLITEM(MAT, VEC, ODE, TIMEDISC, NLTAG) \
    TestCase<MAT, VEC, ODE, NumLib::TIMEDISC, NumLib::NonlinearSolverTag::NLTAG>
#define TCLSEP ,

typedef ::testing::Types<TESTCASESLIST> TestCases;

#undef TESTCASESLIST
#undef TCLSEP
#undef TCLITEM


template<class TestParams>
class NumLibODEIntTyped : public ::testing::Test
{
public:
    using Matrix   = typename TestParams::Matrix;
    using Vector   = typename TestParams::Vector;
    using ODE      = typename TestParams::ODE;
    using TimeDisc = typename TestParams::TimeDisc;

    static const NumLib::NonlinearSolverTag NLTag = TestParams::NLTag;

    static void test()
    {
        for (auto num_timesteps : { 10, 100, 1000 }) {
            run_test_case<Matrix, Vector, TimeDisc, ODE, NLTag>(num_timesteps, TestParams::name);
        }
    }
};


TYPED_TEST_CASE(NumLibODEIntTyped, TestCases);

TYPED_TEST(NumLibODEIntTyped, T1)
{
    TestFixture::test();
}


TEST(NumLibODEInt, ODE3)
{
    const char* name = "dummy";

    // only make sure ODE3 compiles
    run_test_case<EDMatrix, EVector, NumLib::BackwardEuler<EVector>,
                  ODE3<EDMatrix, EVector>,
                  NumLib::NonlinearSolverTag::Newton>
            (0u, name);

    run_test_case<GMatrix, GVector, NumLib::BackwardEuler<GVector>,
                  ODE3<GMatrix, GVector>,
                  NumLib::NonlinearSolverTag::Newton>
            (0u, name);
}



/* TODO Other possible test cases:
 *
 * * check that results are within a specified tolerance
 * * check that results are equal for different matrix/vector types
 * * check that resutls are very close for Picard/Newton
 * * check that the order of time discretization scales correctly
 *   with the timestep size
 */
