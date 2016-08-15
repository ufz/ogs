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

struct Solution
{
    std::vector<double> ts;
    std::vector<GlobalVector> solutions;
};

template<NumLib::NonlinearSolverTag NLTag>
class TestOutput
{
public:
    using TimeDisc = NumLib::TimeDiscretization;
    using NLSolver = NumLib::NonlinearSolver<NLTag>;

    explicit TestOutput(const char* name)
        : _file_name_part(name)
    {}

    template <class ODE>
    Solution run_test(ODE& ode, TimeDisc& timeDisc,
                              const unsigned num_timesteps)
    {
        using ODE_ = ODE;
        using ODET = ODETraits<ODE>;

        Solution sol;

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
        // init_file(delta_t);

        // initial condition
        GlobalVector x0(ode.getMatrixSpecifications().nrows);
        ODET::setIC(x0);

        // write(t0, x0, x0);
        sol.ts.push_back(t0);
        sol.solutions.push_back(x0);

        auto cb = [this, &sol](const double t, GlobalVector const& x) {
            // loopCallback<ODE>(t, x);
            sol.ts.push_back(t);
            sol.solutions.push_back(x);
        };

        if (num_timesteps > 0)
            EXPECT_TRUE(loop.loop(t0, x0, t_end, delta_t, cb));

        return sol;
    }

private:
    /*
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
    */

    const std::string _file_name_part;
    std::unique_ptr<std::ofstream> _file;

    const double _tol = 1e-9;
    const unsigned _maxiter = 20;
};

template <typename TimeDisc, typename ODE, NumLib::NonlinearSolverTag NLTag>
typename std::enable_if<std::is_same<TimeDisc, NumLib::BackwardEuler>::value,
                        Solution>::type
run_test_case(const unsigned num_timesteps, const char* name)
{
    ODE ode;
    TimeDisc timeDisc;

    TestOutput<NLTag> test(name);
    return test.run_test(ode, timeDisc, num_timesteps);
}

template <typename TimeDisc, typename ODE, NumLib::NonlinearSolverTag NLTag>
typename std::enable_if<std::is_same<TimeDisc, NumLib::ForwardEuler>::value,
                        Solution>::type
run_test_case(const unsigned num_timesteps, const char* name)
{
    ODE ode;
    TimeDisc timeDisc;

    TestOutput<NLTag> test(name);
    return test.run_test(ode, timeDisc, num_timesteps);
}

template <typename TimeDisc, typename ODE, NumLib::NonlinearSolverTag NLTag>
typename std::enable_if<std::is_same<TimeDisc, NumLib::CrankNicolson>::value,
                        Solution>::type
run_test_case(const unsigned num_timesteps, const char* name)
{
    ODE ode;
    TimeDisc timeDisc(0.5);

    TestOutput<NLTag> test(name);
    return test.run_test(ode, timeDisc, num_timesteps);
}

template <typename TimeDisc, typename ODE, NumLib::NonlinearSolverTag NLTag>
typename std::enable_if<
    std::is_same<TimeDisc, NumLib::BackwardDifferentiationFormula>::value,
    Solution>::type
run_test_case(const unsigned num_timesteps, const char* name)
{
    ODE ode;
    TimeDisc timeDisc(3);

    TestOutput<NLTag> test(name);
    return test.run_test(ode, timeDisc, num_timesteps);
}

// This class is only here s.t. I don't have to put the members into
// the definition of the macro TCLITEM below.
template <class ODE_, class TimeDisc_>
struct TestCaseBase
{
    using ODE = ODE_;
    using TimeDisc = TimeDisc_;
};

template <class ODE_, class TimeDisc_>
struct TestCase;


// /////////////////////////////////////
//
//  Put new test cases to that list
//
// /////////////////////////////////////
#define TESTCASESLIST \
    TCLITEM(ODE1, BackwardEuler                 , 1e-15  , 0.2) TCLSEP \
    TCLITEM(ODE1, ForwardEuler                  , 1e-13  , 0.22) TCLSEP \
    TCLITEM(ODE1, CrankNicolson                 , 4e-15  , 2.1e-3) TCLSEP \
    TCLITEM(ODE1, BackwardDifferentiationFormula, 4e-15  , 3e-3) TCLSEP \
    \
    TCLITEM(ODE2, BackwardEuler                 , 1.5e-10, 2e-3) TCLSEP \
    TCLITEM(ODE2, ForwardEuler                  , 2e-3   , 2e-3) TCLSEP \
    TCLITEM(ODE2, CrankNicolson                 , 1.5e-10, 8e-6) TCLSEP \
    TCLITEM(ODE2, BackwardDifferentiationFormula, 1.5e-10, 2e-4) TCLSEP \
    \
    TCLITEM(ODE3, BackwardEuler                 , 1e-9   , 0.028) TCLSEP \
    TCLITEM(ODE3, ForwardEuler                  , 1e-13  , 0.03) TCLSEP \
    TCLITEM(ODE3, CrankNicolson                 , 2e-9   , 6e-4) TCLSEP \
    TCLITEM(ODE3, BackwardDifferentiationFormula, 2e-9   , 7e-4)

#define TCLITEM(ODE, TIMEDISC, TOL_PICARD_NEWTON, TOL_ANALYT)                \
    template <>                                                              \
    struct TestCase<ODE, NumLib::TIMEDISC>                                   \
        : TestCaseBase<ODE, NumLib::TIMEDISC> {                              \
        static const char name[];                                            \
        static const double tol_picard_newton;                               \
        static const double tol_analyt;                                      \
    };                                                                       \
    const char TestCase<ODE, NumLib::TIMEDISC>::name[] = #ODE "_" #TIMEDISC; \
    const double TestCase<ODE, NumLib::TIMEDISC>::tol_picard_newton =        \
        (TOL_PICARD_NEWTON);                                                 \
    const double TestCase<ODE, NumLib::TIMEDISC>::tol_analyt = (TOL_ANALYT);
#define TCLSEP

TESTCASESLIST

#undef TCLITEM
#undef TCLSEP

#define TCLITEM(ODE, TIMEDISC, TOL_PICARD_NEWTON, TOL_ANALYT) \
    TestCase<ODE, NumLib::TIMEDISC>
#define TCLSEP ,

typedef ::testing::Types<TESTCASESLIST> TestCases;

#undef TESTCASESLIST
#undef TCLSEP
#undef TCLITEM


template<class TestParams>
class NumLibODEIntTyped : public ::testing::Test
{
public:
    using ODE      = typename TestParams::ODE;
    using TimeDisc = typename TestParams::TimeDisc;

    static void test()
    {
        const unsigned num_timesteps = 100;

        auto const sol_picard =
            run_test_case<TimeDisc, ODE, NumLib::NonlinearSolverTag::Picard>(
                num_timesteps, TestParams::name);

        auto const sol_newton =
            run_test_case<TimeDisc, ODE, NumLib::NonlinearSolverTag::Newton>(
                num_timesteps, TestParams::name);

        // const double tol_picard_newton = 2e-9;

        ASSERT_EQ(sol_picard.ts.size(), sol_newton.ts.size());
        for (std::size_t i = 0; i < sol_picard.ts.size(); ++i) {
            ASSERT_EQ(sol_picard.ts[i], sol_newton.ts[i]);
            for (std::size_t comp = 0; comp < sol_picard.solutions[i].size();
                 ++comp) {
                EXPECT_NEAR(sol_picard.solutions[i][comp],
                            sol_newton.solutions[i][comp],
                            TestParams::tol_picard_newton);

                auto const t = sol_picard.ts[i];
                auto const sol_analyt = ODETraits<ODE>::solution(t);

                EXPECT_NEAR(sol_picard.solutions[i][comp], sol_analyt[comp],
                            TestParams::tol_analyt);
                EXPECT_NEAR(sol_newton.solutions[i][comp], sol_analyt[comp],
                            TestParams::tol_analyt);
            }
        }
    }
};

TYPED_TEST_CASE(NumLibODEIntTyped, TestCases);

TYPED_TEST(NumLibODEIntTyped, T1)
{
    TestFixture::test();
}


/* TODO Other possible test cases:
 *
 * * check that results are within a specified tolerance
 * * check that results are equal for different matrix/vector types
 * * check that resutls are very close for Picard/Newton
 * * check that the order of time discretization scales correctly
 *   with the timestep size
 */
