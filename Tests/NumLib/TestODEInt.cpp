/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gtest/gtest.h>

#include <fstream>
#include <memory>
#include <typeinfo>

#include <logog/include/logog.hpp>

#include "BaseLib/BuildInfo.h"
#include "NumLib/NumericsConfig.h"
#include "NumLib/ODESolver/TimeLoopSingleODE.h"
#include "NumLib/ODESolver/ConvergenceCriterionDeltaX.h"
#include "ODEs.h"

using GMatrix = GlobalMatrix;
using GVector = GlobalVector;


template<typename Matrix, typename Vector, NumLib::NonlinearSolverTag NLTag>
class TestOutput
{
public:
    using TimeDisc = NumLib::TimeDiscretization;
    using NLSolver = NumLib::NonlinearSolver<NLTag>;

    explicit TestOutput(const char* name)
        : _file_name_part(name)
    {}

    template<typename ODE>
    void run_test(ODE& ode, TimeDisc& timeDisc)
    {
        run_test<ODE>(ode, timeDisc, 10); // by default make 10 timesteps
    }

    template <class ODE>
    void run_test(ODE& ode, TimeDisc& timeDisc, const unsigned num_timesteps)
    {
        using ODE_ = ODE;
        using ODET = ODETraits<ODE>;

        NumLib::TimeDiscretizedODESystem<ODE_::ODETag, NLTag>
                ode_sys(ode, timeDisc);

        auto linear_solver = std::unique_ptr<GlobalLinearSolver>{
            new GlobalLinearSolver{"", nullptr}};
        auto conv_crit = std::unique_ptr<NumLib::ConvergenceCriterion>(
            new NumLib::ConvergenceCriterionDeltaX(
                _tol, boost::none, MathLib::VecNormType::NORM2));
        std::unique_ptr<NLSolver> nonlinear_solver(
            new NLSolver(*linear_solver, _maxiter));

        NumLib::TimeLoopSingleODE<NLTag> loop(ode_sys, std::move(linear_solver),
                                              std::move(nonlinear_solver),
                                              std::move(conv_crit));

        const double t0      = ODET::t0;
        const double t_end   = ODET::t_end;
        const double delta_t = (num_timesteps == 0) ? -1.0
                                                    : ((t_end-t0) / num_timesteps);

        INFO("Running test %s with %u timesteps of size %g s.",
             _file_name_part.c_str(), num_timesteps, delta_t);
        init_file(delta_t);

        // initial condition
        Vector x0(ode.getMatrixSpecifications().nrows);
        ODET::setIC(x0);

        write(t0, x0, x0);

        auto cb = [this](const double t, Vector const& x) {
            loopCallback<ODE>(t, x);
        };

        if (num_timesteps > 0)
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
        for (decltype(x_num.size()) i = 0; i < x_num.size(); ++i)
            *_file << '\t' << x_num[i];
        for (decltype(x_ana.size()) i = 0; i < x_ana.size(); ++i)
            *_file << '\t' << x_ana[i];
        *_file << "\n";
    }

    template <class Ode>
    void loopCallback(const double t, Vector const& x)
    {
        write(t, x, ODETraits<Ode>::solution(t));
    }

    const std::string _file_name_part;
    std::unique_ptr<std::ofstream> _file;

    const double _tol = 1e-9;
    const unsigned _maxiter = 20;
};


template<typename Matrix, typename Vector, typename TimeDisc, typename ODE,
         NumLib::NonlinearSolverTag NLTag>
typename std::enable_if<std::is_same<TimeDisc, NumLib::BackwardEuler >::value>::type
run_test_case(const unsigned num_timesteps, const char* name)
{
    ODE ode;
    TimeDisc timeDisc;

    TestOutput<Matrix, Vector, NLTag> test(name);
    test.run_test(ode, timeDisc, num_timesteps);
}

template<typename Matrix, typename Vector, typename TimeDisc, typename ODE,
         NumLib::NonlinearSolverTag NLTag>
typename std::enable_if<std::is_same<TimeDisc, NumLib::ForwardEuler >::value>::type
run_test_case(const unsigned num_timesteps, const char* name)
{
    ODE ode;
    TimeDisc timeDisc;

    TestOutput<Matrix, Vector, NLTag> test(name);
    test.run_test(ode, timeDisc, num_timesteps);
}

template<typename Matrix, typename Vector, typename TimeDisc, typename ODE,
         NumLib::NonlinearSolverTag NLTag>
typename std::enable_if<std::is_same<TimeDisc, NumLib::CrankNicolson >::value>::type
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
    std::is_same<TimeDisc, NumLib::BackwardDifferentiationFormula >::value>::type
run_test_case(const unsigned num_timesteps, const char* name)
{
    ODE ode;
    TimeDisc timeDisc(3);

    TestOutput<Matrix, Vector, NLTag> test(name);
    test.run_test(ode, timeDisc, num_timesteps);
}


// This class is only here s.t. I don't have to put the members into
// the definition of the macro TCLITEM below.
template <typename Matrix_, typename Vector_, class ODE_, class TimeDisc_,
          NumLib::NonlinearSolverTag NLTag_>
struct TestCaseBase
{
    using Matrix = Matrix_;
    using Vector = Vector_;
    using ODE = ODE_;
    using TimeDisc = TimeDisc_;
    static const NumLib::NonlinearSolverTag NLTag = NLTag_;
};

template <typename Matrix_, typename Vector_, class ODE_, class TimeDisc_,
          NumLib::NonlinearSolverTag NLTag_>
struct TestCase;


// /////////////////////////////////////
//
//  Put new test cases to that list
//
// /////////////////////////////////////
#define TESTCASESLIST \
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
        for (auto num_timesteps : { 10u, 100u })
        {
            run_test_case<Matrix, Vector, TimeDisc, ODE, NLTag>(
                        num_timesteps, TestParams::name);
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
    {
        run_test_case<GMatrix, GVector, NumLib::BackwardEuler, ODE3,
                      NumLib::NonlinearSolverTag::Newton>(0u, name);
    }
}



/* TODO Other possible test cases:
 *
 * * check that results are within a specified tolerance
 * * check that results are equal for different matrix/vector types
 * * check that resutls are very close for Picard/Newton
 * * check that the order of time discretization scales correctly
 *   with the timestep size
 */
