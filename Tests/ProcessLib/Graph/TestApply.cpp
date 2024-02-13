/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gtest/gtest.h>

#include "BaseLib/StrongType.h"
#include "ProcessLib/Graph/Apply.h"

namespace
{
void f_zero_args() {}
}  // namespace

namespace PG = ProcessLib::Graph;
namespace PGD = ProcessLib::Graph::detail;

class ProcessLibGraphApply : public testing::Test
{
protected:
    using SInt = BaseLib::StrongType<int, struct SIntTag>;
    using SDouble = BaseLib::StrongType<double, struct SDoubleTag>;

    // --- tested functions ----------------------------------------------------

    // --- single output argument

    static constexpr auto lambda_increment = [](int& i) { ++i; };

    struct SIncrementF
    {
        void f(int& i) const { ++i; }
    };
    struct SIncrementEval
    {
        void eval(int& i) const { ++i; }
    };

    // returns a vector of functions running the code under test
    template <typename Data>
    static std::vector<std::function<void(Data&)>>
    getRunnersSIncrementSingleTuple()
    {
        return {[](Data& args) { PGD::applyImpl(lambda_increment, args.t); },
                // ---
                [](Data& args)
                {
                    SIncrementF const s;
                    PGD::applyImpl(&SIncrementF::f, s, args.t);
                },
                // ---
                [](Data& args) { PG::apply(lambda_increment, args.t); },
                // ---
                [](Data& args)
                {
                    SIncrementEval const s;
                    PG::eval(s, args.t);
                }};
    }

    // --- single argument and return value

    struct SIdentity
    {
        double operator()(double d) { return d; }
        double eval(double d) { return d; }
        double f(double d) { return d; }
    };

    // returns a vector of functions running the code under test
    template <typename Data>
    static std::vector<std::function<double(Data&)>>
    getRunnersSIdentitySingleTuple()
    {
        return {[](Data& args) { return PGD::applyImpl(SIdentity{}, args.t); },
                // ---
                [](Data& args)
                {
                    SIdentity s;
                    return PGD::applyImpl(&SIdentity::f, s, args.t);
                },
                // ---
                [](Data& args)
                {
                    SIdentity s;
                    return PG::apply(s, args.t);
                },
                // ---
                [](Data& args)
                {
                    SIdentity s;
                    return PG::eval(s, args.t);
                }};
    }

    // --- input and output argument

    struct SAddSecondToFirstF
    {
        void f(double& d, SInt const& si) const { d += *si; }
    };
    struct SAddSecondToFirstFunctionObject
    {
        void operator()(double& d, SInt const& si) const { d += *si; }
    };
    struct SAddSecondToFirstEval
    {
        void eval(double& d, SInt const& si) const { d += *si; }
    };

    // returns a vector of functions running the code under test
    template <typename Data>
    static std::vector<std::function<void(Data&)>>
    getRunnersSAddSecondToFirstSingleTuple()
    {
        return {[](Data& args)
                {
                    SAddSecondToFirstF const s;
                    PGD::applyImpl(&SAddSecondToFirstF::f, s, args.t);
                },
                // ---
                [](Data& args)
                {
                    SAddSecondToFirstFunctionObject const s;
                    PGD::applyImpl(s, args.t);
                },
                // ---
                [](Data& args)
                {
                    SAddSecondToFirstFunctionObject const s;
                    PG::apply(s, args.t);
                },
                // ---
                [](Data& args)
                {
                    SAddSecondToFirstEval const s;
                    PG::eval(s, args.t);
                }};
    }

    // --- multiple arguments

    struct SMultipleArguments
    {
        int f(SInt const& si, double const& d, SDouble& sd, int i)
        {
            *sd += *si + d - i;
            return 2 * i;
        }
        int operator()(SInt const& si, double const& d, SDouble& sd, int i)
        {
            return f(si, d, sd, i);
        }
        int eval(SInt const& si, double const& d, SDouble& sd, int i)
        {
            return f(si, d, sd, i);
        }
    };

    // returns a vector of functions running the code under test
    template <typename Data>
    static std::vector<std::function<int(Data&)>>
    getRunnersSMultipleArgumentsSingleTuple()
    {
        return {[](Data& args)
                {
                    SMultipleArguments s;
                    return PGD::applyImpl(&SMultipleArguments::f, s, args.t);
                },
                // ---
                [](Data& args)
                { return PGD::applyImpl(SMultipleArguments{}, args.t); },
                // ---
                [](Data& args)
                {
                    SMultipleArguments s;
                    return PG::apply(s, args.t);
                },
                // ---
                [](Data& args)
                {
                    SMultipleArguments s;
                    return PG::eval(s, args.t);
                }};
    }

    // --- test data -----------------------------------------------------------

    struct DataSingleTuple
    {
        std::tuple<int, double, SInt const, SDouble> t{5, 6.5, 2, 0.25};

        // keeps original data to compute reference results for checks
        std::tuple<int, double, SInt, SDouble> const init{t};
    };

    struct DataSingleTupleWithReferences
    {
        double d = 6.5;
        SInt si{2};
        SDouble sd{0.25};

        std::tuple<int, double&, SInt const&, SDouble&> t{5, d, si, sd};

        // keeps original data to compute reference results for checks
        std::tuple<int, double, SInt, SDouble> const init{t};
    };

    struct DataMultipleTuples
    {
        double d = 6.5;
        SInt si{2};
        SDouble sd{0.25};

        std::tuple<SDouble&> t1{sd};
        std::tuple<SInt const&> const t2{si};
        std::tuple<int, double&> t3{5, d};
        std::tuple<> t4;
        std::tuple<std::string> t5{"unrelated"};

        // keeps original data to compute reference results for checks
        std::tuple<SDouble, SInt, int, double, std::string> const init{
            std::tuple_cat(t1, t2, t3, t4, t5)};
    };
};

TEST_F(ProcessLibGraphApply, ZeroArgFct)
{
    struct S
    {
        void operator()() {}
        void f() {}
        int g() { return ++i; }

    private:
        int i = 7;
    };

    struct S2
    {
        int operator()() { return --i; }

    private:
        int i = 20;
    };

    struct SConst
    {
        void operator()() const {}
        void f() const {}
        int g() const { return 41; }
    };

    auto this_should_compile = []<typename... T>(T&&... args)
    {
        std::tuple<> nothing;
        std::tuple<char, int, double> something{'a', 0, 0.5};

        auto const& cnothing = nothing;
        auto const& csomething = something;

        // We test an implementation detail, because it serves as a "backend" to
        // both apply() and eval(), but is more generic and easier to test.
        PGD::applyImpl(std::forward<T>(args)...);
        PGD::applyImpl(std::forward<T>(args)..., nothing);
        PGD::applyImpl(std::forward<T>(args)..., cnothing);
        PGD::applyImpl(std::forward<T>(args)..., something);
        PGD::applyImpl(std::forward<T>(args)..., csomething);
        return PGD::applyImpl(std::forward<T>(args)..., nothing, csomething,
                              cnothing);
    };

    {
        auto lambda = []() {};
        auto lambda_mut = []() mutable {};
        auto lambda_capture = [i = 5]() { return i; };
        auto lambda_capture_mut = [i = 13]() mutable { return ++i; };

        this_should_compile(lambda);                        // lambda
        this_should_compile(lambda_mut);                    // mutable ...
        EXPECT_EQ(5, this_should_compile(lambda_capture));  // ... with capture

        lambda_capture_mut();  // mutable ... with capture
        EXPECT_EQ(13 + 6 + 1, this_should_compile(lambda_capture_mut))
            << "lambda_capture_mut was expected to be called six times and "
               "increase its internal counter. This test guarantees that the "
               "implementation does not copy the function object.";
    }

    this_should_compile(f_zero_args);  // free function

    {
        S s;
        this_should_compile(s);         // function object (void)
        this_should_compile(&S::f, s);  // void member function

        s.g();  // non-void member function
        EXPECT_EQ(7 + 6 + 1, this_should_compile(&S::g, s))
            << "S::g() was expected to be called six times and increase its "
               "internal counter. This test guarantees that the implementation "
               "does not copy the object.";
    }

    {
        S2 s2;  // function object (non-void)
        s2();
        EXPECT_EQ(20 - 6 - 1, this_should_compile(s2))
            << "S2::operator()() was expected to be called six times and "
               "decrease its internal counter. This test guarantees that the "
               "implementation does not copy the function object.";
    }

    {
        SConst const sconst;
        this_should_compile(sconst);  // const function object
        this_should_compile(&SConst::f,
                            sconst);  // const member function (void)
        EXPECT_EQ(41,
                  this_should_compile(&SConst::g, sconst));  // ... (non-void)
    }

    // temporary lambda
    EXPECT_EQ(-8 - 6, this_should_compile([i = -8]() mutable { return --i; }))
        << "The lambda was expected to be called six times and decrease "
           "its internal counter. This test guarantees that the implementation "
           "does not copy the function object.";

    // temporary function objects
    this_should_compile(S{});
    EXPECT_EQ(20 - 6, this_should_compile(S2{}));
    this_should_compile(SConst{});

    // member functions and temporary objects
    this_should_compile(&S::f, S{});
    EXPECT_EQ(7 + 6, this_should_compile(&S::g, S{}));
    this_should_compile(&SConst::f, SConst{});
    EXPECT_EQ(41, this_should_compile(&SConst::g, SConst{}));
}

TEST_F(ProcessLibGraphApply, SingleTupleSingleOutputArgument)
{
    auto check = [](DataSingleTuple const& args)
    {
        EXPECT_EQ(get<int>(args.init) + 1, get<int>(args.t));
        EXPECT_EQ(get<double>(args.init), get<double>(args.t));
        EXPECT_EQ(*get<SInt>(args.init), *get<SInt const>(args.t));
        EXPECT_EQ(*get<SDouble>(args.init), *get<SDouble>(args.t));
    };

    for (auto& runner : getRunnersSIncrementSingleTuple<DataSingleTuple>())
    {
        DataSingleTuple args;
        runner(args);
        check(args);
    }
}

TEST_F(ProcessLibGraphApply, SingleTupleSingleArgumentAndReturnValue)
{
    auto check = [](DataSingleTuple const& args)
    {
        EXPECT_EQ(get<int>(args.init), get<int>(args.t));
        EXPECT_EQ(get<double>(args.init), get<double>(args.t));
        EXPECT_EQ(*get<SInt>(args.init), *get<SInt const>(args.t));
        EXPECT_EQ(*get<SDouble>(args.init), *get<SDouble>(args.t));
    };

    for (auto& runner : getRunnersSIdentitySingleTuple<DataSingleTuple>())
    {
        DataSingleTuple args;

        EXPECT_EQ(get<double>(args.init), runner(args));

        check(args);
    }
}

TEST_F(ProcessLibGraphApply, SingleTupleSingleArgumentAndReturnReference)
{
    auto check = [](DataSingleTuple const& args)
    {
        EXPECT_EQ(get<int>(args.init), get<int>(args.t));
        EXPECT_EQ(get<double>(args.init), get<double>(args.t));
        EXPECT_EQ(*get<SInt>(args.init), *get<SInt const>(args.t));
        EXPECT_EQ(*get<SDouble>(args.init), *get<SDouble>(args.t));
    };

    auto lambda = [](double& d) -> double& { return d; };
    auto lambda_capture = [i = 1.0](double& d) -> double&
    {
        d *= i;
        return d;
    };
    auto lambda_mut = [](double& d) mutable -> double& { return d; };
    auto lambda_capture_mut = [i = 1.0](double& d) mutable -> double&
    {
        d *= i;
        return d;
    };

    struct S
    {
        double& f(double& d) { return d; }
        double& operator()(double& d) { return d; }
    };

    struct SConst
    {
        double& f(double& d) const { return d; }
        double& operator()(double& d) const { return d; }
    };

    std::function<double&(DataSingleTuple&)> runners[] = {
        [&lambda](DataSingleTuple& args) -> double&
        { return PGD::applyImpl(lambda, args.t); },
        // ---
        [&lambda_capture](DataSingleTuple& args) -> double&
        { return PGD::applyImpl(lambda_capture, args.t); },
        // ---
        [&lambda_mut](DataSingleTuple& args) -> double&
        { return PGD::applyImpl(lambda_mut, args.t); },
        // ---
        [&lambda_capture_mut](DataSingleTuple& args) -> double&
        { return PGD::applyImpl(lambda_capture_mut, args.t); },
        // ---
        [](DataSingleTuple& args) -> double&
        {
            S s;
            return PGD::applyImpl(s, args.t);
        },
        // ---
        [](DataSingleTuple& args) -> double&
        {
            S s;
            return PGD::applyImpl(&S::f, s, args.t);
        },
        // ---
        [](DataSingleTuple& args) -> double&
        {
            SConst s;
            return PGD::applyImpl(s, args.t);
        },
        // ---
        [](DataSingleTuple& args) -> double&
        {
            SConst s;
            return PGD::applyImpl(&SConst::f, s, args.t);
        }};

    for (auto& runner : runners)
    {
        DataSingleTuple args;

        auto& res = runner(args);
        EXPECT_EQ(get<double>(args.init), res);

        check(args);

        // check that the returned reference points to the right data
        EXPECT_EQ(&get<double>(args.t), &res);
        EXPECT_NE(31, res)
            << "We want to modify res to that value in the next step.";
        res = 31;
        EXPECT_EQ(31, get<double>(args.t));
    }
}

TEST_F(ProcessLibGraphApply, SingleTupleInputAndOutputArgument)
{
    auto check = [](DataSingleTuple const& args)
    {
        EXPECT_EQ(get<int>(args.init), get<int>(args.t));
        EXPECT_EQ(get<double>(args.init) + *get<SInt>(args.init),
                  get<double>(args.t));
        EXPECT_EQ(*get<SInt>(args.init), *get<SInt const>(args.t));
        EXPECT_EQ(*get<SDouble>(args.init), *get<SDouble>(args.t));
    };

    for (auto& runner :
         getRunnersSAddSecondToFirstSingleTuple<DataSingleTuple>())
    {
        DataSingleTuple args;
        runner(args);
        check(args);
    }
}

TEST_F(ProcessLibGraphApply, SingleTupleMultipleArguments)
{
    auto check = [](DataSingleTuple const& args)
    {
        EXPECT_EQ(get<int>(args.init), get<int>(args.t));
        EXPECT_EQ(get<double>(args.init), get<double>(args.t));
        EXPECT_EQ(*get<SInt>(args.init), *get<SInt const>(args.t));
        EXPECT_EQ(*get<SDouble>(args.init) - get<int>(args.init) +
                      get<double>(args.init) + *get<SInt>(args.init),
                  *get<SDouble>(args.t));
    };

    for (auto& runner :
         getRunnersSMultipleArgumentsSingleTuple<DataSingleTuple>())
    {
        DataSingleTuple args;

        EXPECT_EQ(2 * get<int>(args.init), runner(args));

        check(args);
    }
}

TEST_F(ProcessLibGraphApply, SingleTupleWithReferencesSingleOutputArgument)
{
    auto check = [](DataSingleTupleWithReferences const& args)
    {
        EXPECT_EQ(get<int>(args.init) + 1, get<int>(args.t));
        EXPECT_EQ(get<double>(args.init), get<double&>(args.t));
        EXPECT_EQ(*get<SInt>(args.init), *get<SInt const&>(args.t));
        EXPECT_EQ(*get<SDouble>(args.init), *get<SDouble&>(args.t));
    };

    for (auto& runner :
         getRunnersSIncrementSingleTuple<DataSingleTupleWithReferences>())
    {
        DataSingleTupleWithReferences args;
        runner(args);
        check(args);
    }
}

TEST_F(ProcessLibGraphApply,
       SingleTupleWithReferencesSingleArgumentAndReturnValue)
{
    auto check = [](DataSingleTupleWithReferences const& args)
    {
        EXPECT_EQ(get<int>(args.init), get<int>(args.t));
        EXPECT_EQ(get<double>(args.init), get<double&>(args.t));
        EXPECT_EQ(*get<SInt>(args.init), *get<SInt const&>(args.t));
        EXPECT_EQ(*get<SDouble>(args.init), *get<SDouble&>(args.t));
    };

    for (auto& runner :
         getRunnersSIdentitySingleTuple<DataSingleTupleWithReferences>())
    {
        DataSingleTupleWithReferences args;

        EXPECT_EQ(get<double>(args.init), runner(args));

        check(args);
    }
}

TEST_F(ProcessLibGraphApply, SingleTupleWithReferencesInputAndOutputArgument)
{
    auto check = [](DataSingleTupleWithReferences const& args)
    {
        EXPECT_EQ(get<int>(args.init), get<int>(args.t));
        EXPECT_EQ(get<double>(args.init) + *get<SInt>(args.init),
                  get<double&>(args.t));
        EXPECT_EQ(get<double>(args.init) + *get<SInt>(args.init), args.d)
            << "The original variable has not been modified correctly.";
        EXPECT_EQ(*get<SInt>(args.init), *get<SInt const&>(args.t));
        EXPECT_EQ(*get<SDouble>(args.init), *get<SDouble&>(args.t));
    };

    for (auto& runner : getRunnersSAddSecondToFirstSingleTuple<
             DataSingleTupleWithReferences>())
    {
        DataSingleTupleWithReferences args;
        runner(args);
        check(args);
    }
}

TEST_F(ProcessLibGraphApply, SingleTupleWithReferencesMultipleArguments)
{
    auto check = [](DataSingleTupleWithReferences const& args)
    {
        EXPECT_EQ(get<int>(args.init), get<int>(args.t));
        EXPECT_EQ(get<double>(args.init), get<double&>(args.t));
        EXPECT_EQ(*get<SInt>(args.init), *get<SInt const&>(args.t));
        EXPECT_EQ(*get<SDouble>(args.init) - get<int>(args.init) +
                      get<double>(args.init) + *get<SInt>(args.init),
                  *get<SDouble&>(args.t));
        EXPECT_EQ(*get<SDouble>(args.init) - get<int>(args.init) +
                      get<double>(args.init) + *get<SInt>(args.init),
                  *args.sd)
            << "The original variable has not been modified correctly.";
    };

    for (auto& runner : getRunnersSMultipleArgumentsSingleTuple<
             DataSingleTupleWithReferences>())
    {
        DataSingleTupleWithReferences args;

        EXPECT_EQ(2 * get<int>(args.init), runner(args));

        check(args);
    }
}

TEST_F(ProcessLibGraphApply, MultipleTuplesSingleOutputArgument)
{
    auto check = [](DataMultipleTuples const& args)
    {
        EXPECT_EQ(get<int>(args.init) + 1, get<int>(args.t3));
        EXPECT_EQ(get<double>(args.init), get<double&>(args.t3));
        EXPECT_EQ(*get<SInt>(args.init), *get<SInt const&>(args.t2));
        EXPECT_EQ(*get<SDouble>(args.init), *get<SDouble&>(args.t1));
        EXPECT_EQ(get<std::string>(args.init), get<std::string>(args.t5));
    };

    std::function<void(DataMultipleTuples&)> runners[] = {
        [](DataMultipleTuples& args)
        {
            PGD::applyImpl(
                lambda_increment, args.t1, args.t2, args.t3, args.t4, args.t5);
        },
        // ---
        [](DataMultipleTuples& args)
        {
            SIncrementF const s;
            PGD::applyImpl(&SIncrementF::f,
                           s,
                           args.t1,
                           args.t2,
                           args.t3,
                           args.t4,
                           args.t5);
        },
        // ---
        [](DataMultipleTuples& args) {
            PG::apply(
                lambda_increment, args.t1, args.t2, args.t3, args.t4, args.t5);
        },
        // ---
        [](DataMultipleTuples& args)
        {
            SIncrementEval const s;
            PG::eval(s, args.t1, args.t2, args.t3, args.t4, args.t5);
        }};

    for (auto& runner : runners)
    {
        DataMultipleTuples args;
        runner(args);
        check(args);
    }
}

TEST_F(ProcessLibGraphApply, MultipleTuplesSingleArgumentAndReturnValue)
{
    auto check = [](DataMultipleTuples const& args)
    {
        EXPECT_EQ(get<int>(args.init), get<int>(args.t3));
        EXPECT_EQ(get<double>(args.init), get<double&>(args.t3));
        EXPECT_EQ(*get<SInt>(args.init), *get<SInt const&>(args.t2));
        EXPECT_EQ(*get<SDouble>(args.init), *get<SDouble&>(args.t1));
        EXPECT_EQ(get<std::string>(args.init), get<std::string>(args.t5));
    };

    std::function<double(DataMultipleTuples&)> runners[] = {
        [](DataMultipleTuples& args)
        {
            return PGD::applyImpl(
                SIdentity{}, args.t2, args.t3, args.t4, args.t1, args.t5);
        },
        // ---
        [](DataMultipleTuples& args)
        {
            SIdentity s;
            return PGD::applyImpl(
                &SIdentity::f, s, args.t2, args.t3, args.t4, args.t1, args.t5);
        },
        // ---
        [](DataMultipleTuples& args)
        {
            SIdentity s;
            return PG::apply(s, args.t2, args.t3, args.t4, args.t1, args.t5);
        },
        // ---
        [](DataMultipleTuples& args)
        {
            SIdentity s;
            return PG::eval(s, args.t2, args.t3, args.t4, args.t1, args.t5);
        }};

    for (auto& runner : runners)
    {
        DataMultipleTuples args;

        EXPECT_EQ(get<double>(args.init), runner(args));

        check(args);
    }
}

TEST_F(ProcessLibGraphApply, MultipleTuplesInputAndOutputArgument)
{
    auto check = [](DataMultipleTuples const& args)
    {
        EXPECT_EQ(get<int>(args.init), get<int>(args.t3));
        EXPECT_EQ(get<double>(args.init) + *get<SInt>(args.init),
                  get<double&>(args.t3));
        EXPECT_EQ(get<double>(args.init) + *get<SInt>(args.init), args.d)
            << "The original variable has not been modified correctly.";
        EXPECT_EQ(*get<SInt>(args.init), *get<SInt const&>(args.t2));
        EXPECT_EQ(*get<SDouble>(args.init), *get<SDouble&>(args.t1));
        EXPECT_EQ(get<std::string>(args.init), get<std::string>(args.t5));
    };

    std::function<void(DataMultipleTuples&)> runners[] = {
        [](DataMultipleTuples& args)
        {
            SAddSecondToFirstF const s;
            PGD::applyImpl(&SAddSecondToFirstF::f,
                           s,
                           args.t5,
                           args.t4,
                           args.t3,
                           args.t2,
                           args.t1);
        },
        // ---
        [](DataMultipleTuples& args)
        {
            SAddSecondToFirstFunctionObject const s;
            PGD::applyImpl(s, args.t5, args.t4, args.t3, args.t2, args.t1);
        },
        // ---
        [](DataMultipleTuples& args)
        {
            SAddSecondToFirstFunctionObject const s;
            PG::apply(s, args.t5, args.t4, args.t3, args.t2, args.t1);
        },
        // ---
        [](DataMultipleTuples& args)
        {
            SAddSecondToFirstEval const s;
            PG::eval(s, args.t5, args.t4, args.t3, args.t2, args.t1);
        }};

    for (auto& runner : runners)
    {
        DataMultipleTuples args;
        runner(args);
        check(args);
    }
}

TEST_F(ProcessLibGraphApply, MultipleTuplesMultipleArguments)
{
    auto check = [](DataMultipleTuples const& args)
    {
        EXPECT_EQ(get<int>(args.init), get<int>(args.t3));
        EXPECT_EQ(get<double>(args.init), get<double&>(args.t3));
        EXPECT_EQ(*get<SInt>(args.init), *get<SInt const&>(args.t2));
        EXPECT_EQ(*get<SDouble>(args.init) - get<int>(args.init) +
                      get<double>(args.init) + *get<SInt>(args.init),
                  *get<SDouble&>(args.t1));
        EXPECT_EQ(*get<SDouble>(args.init) - get<int>(args.init) +
                      get<double>(args.init) + *get<SInt>(args.init),
                  *args.sd)
            << "The original variable has not been modified correctly.";
        EXPECT_EQ(get<std::string>(args.init), get<std::string>(args.t5));
    };

    std::function<int(DataMultipleTuples&)> runners[] = {
        [](DataMultipleTuples& args)
        {
            SMultipleArguments s;
            return PGD::applyImpl(&SMultipleArguments::f,
                                  s,
                                  args.t1,
                                  args.t5,
                                  args.t2,
                                  args.t3,
                                  args.t4);
        },
        // ---
        [](DataMultipleTuples& args)
        {
            return PGD::applyImpl(SMultipleArguments{},
                                  args.t1,
                                  args.t5,
                                  args.t2,
                                  args.t3,
                                  args.t4);
        },
        // ---
        [](DataMultipleTuples& args)
        {
            SMultipleArguments s;
            return PG::apply(s, args.t1, args.t5, args.t2, args.t3, args.t4);
        },
        // ---
        [](DataMultipleTuples& args)
        {
            SMultipleArguments s;
            return PG::eval(s, args.t1, args.t5, args.t2, args.t3, args.t4);
        }};

    for (auto& runner : runners)
    {
        DataMultipleTuples args;

        EXPECT_EQ(2 * get<int>(args.init), runner(args));

        check(args);
    }
}
