/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gtest/gtest.h>
#include <logog/include/logog.hpp>

#include "BaseLib/Functional.h"
#include "NumLib/NamedFunctionProvider.h"
#include "NumLib/NamedFunctionCaller.h"
#include "Tests/InstanceCounter.h"

class F : public NumLib::NamedFunctionProvider
{
public:
    double f(double arg_g, double arg_y) const
    {
        return arg_g + arg_y;
    }

    std::vector<NumLib::NamedFunction>
    getNamedFunctions() const override
    {
        return {{"f", {"g_arg", "y"}, BaseLib::easyBind(&F::f, this)}};
    }
};

class G : public NumLib::NamedFunctionProvider
{
public:
    double g(double arg_x) const
    {
        return -arg_x;
    }

    std::vector<NumLib::NamedFunction>
    getNamedFunctions() const override
    {
        return {{"g", {"x"}, BaseLib::easyBind(&G::g, this)}};
    }
};

class H : public InstanceCounter<H>
{
public:
    H(double const z) : _z(z) {}
    double h(double arg_x, double arg_y) { return arg_x * arg_y - _z; }
    double setZ(double const z)
    {
        _z = z;
        return z;
    }

private:
    double _z;
};

class I : public NumLib::NamedFunctionProvider
{
public:
    double i(double arg_x) const
    {
        return -arg_x;
    }

    std::vector<NumLib::NamedFunction>
    getNamedFunctions() const override
    {
        return {{"i", {"x"}, BaseLib::easyBind(&I::i, this)}};
    }
};

TEST(NumLib, NamedFunctionCaller)
{
    F f_inst;
    G g_inst;

    NumLib::NamedFunctionCaller caller{ "x", "y" };

    for (auto&& f_named : f_inst.getNamedFunctions()) {
        caller.addNamedFunction(std::move(f_named));
    }

    caller.plug("f", "g_arg", "g");
    caller.plug("f", "y", "y");
    caller.plug("g", "x", "x");

    // test if adding function after plug works
    for (auto&& g_named : g_inst.getNamedFunctions()) {
        caller.addNamedFunction(std::move(g_named));
    }

    caller.applyPlugs();

    double x = 1.0;
    double y = 2.0;

    auto const g_caller = caller.getSpecificFunctionCaller("g");
    DBUG("calling %s", caller.getCallExpression("g").c_str());
    EXPECT_EQ(g_inst.g(x), g_caller.call({x, y}));

    auto const f_caller = caller.getSpecificFunctionCaller("f");
    DBUG("calling %s", caller.getCallExpression("f").c_str());
    EXPECT_EQ(f_inst.f(g_inst.g(x), y), f_caller.call({x, y}));
}

TEST(NumLib, NamedFunctionCallerCyclicGraph)
{
    // Construct a cyclic case with f(g(i(f(...), y)))
    F f_inst;
    G g_inst;
    I i_inst;

    NumLib::NamedFunctionCaller caller{ "x", "y" };

    for (auto&& f_named : f_inst.getNamedFunctions()) {
        caller.addNamedFunction(std::move(f_named));
    }
    for (auto&& g_named : g_inst.getNamedFunctions()) {
        caller.addNamedFunction(std::move(g_named));
    }
    for (auto&& i_named : i_inst.getNamedFunctions()) {
        caller.addNamedFunction(std::move(i_named));
    }

    caller.plug("f", "g_arg", "g");
    caller.plug("f", "y", "y");
    caller.plug("g", "x", "i");
    caller.plug("i", "x", "f");

    ASSERT_ANY_THROW(caller.applyPlugs());
}

TEST(NumLib, NamedFunctionNoLeaks)
{
    auto num_const = InstanceCounter<H>::getNumberOfConstructions();
    auto num_move = InstanceCounter<H>::getNumberOfMoves();
    auto num_copy = InstanceCounter<H>::getNumberOfCopies();
    auto num_inst = InstanceCounter<H>::getNumberOfInstances();

    {
        H h_inst(1.0);
        EXPECT_EQ(num_const+1, InstanceCounter<H>::getNumberOfConstructions());
        EXPECT_EQ(num_inst+1, InstanceCounter<H>::getNumberOfInstances());
        InstanceCounter<H>::update(num_const, num_move, num_copy, num_inst);

        auto h_fct = NumLib::NamedFunction("h", {"x", "y"},
                                           BaseLib::easyBind(&H::h, h_inst));
        EXPECT_EQ(num_const, InstanceCounter<H>::getNumberOfConstructions());
        EXPECT_EQ(num_inst, InstanceCounter<H>::getNumberOfInstances());

        EXPECT_EQ(5.0, h_fct.call({ 2.0, 3.0 }));
        h_inst.setZ(2.0);
        EXPECT_EQ(4.0, h_fct.call({ 2.0, 3.0 }));

        auto h_bind = BaseLib::easyBind(&H::h, H{3.0});
        EXPECT_EQ(num_const+1, InstanceCounter<H>::getNumberOfConstructions());
        EXPECT_EQ(num_inst+1, InstanceCounter<H>::getNumberOfInstances());
        InstanceCounter<H>::update(num_const, num_move, num_copy, num_inst);

        // Move an object. NamedFunction will implicitly do memory management.
        auto h_fct2 = NumLib::NamedFunction("h", {"x", "y"}, std::move(h_bind));
        EXPECT_EQ(num_copy, InstanceCounter<H>::getNumberOfCopies());
        EXPECT_EQ(num_const, InstanceCounter<H>::getNumberOfConstructions());
        EXPECT_EQ(num_inst, InstanceCounter<H>::getNumberOfInstances());

        EXPECT_EQ(5.0, h_fct2.call({ 2.0, 4.0 }));

        // copy
        NumLib::NamedFunction h_fct3 = h_fct2;
        EXPECT_EQ(num_const, InstanceCounter<H>::getNumberOfConstructions());
        EXPECT_EQ(num_copy+1, InstanceCounter<H>::getNumberOfCopies());
        EXPECT_EQ(num_inst+1, InstanceCounter<H>::getNumberOfInstances());
        InstanceCounter<H>::update(num_const, num_move, num_copy, num_inst);

        EXPECT_EQ(-1.0, h_fct3.call({ 2.0, 1.0 }));

        // move
        NumLib::NamedFunction h_fct4 = std::move(h_fct3);
        EXPECT_INSTANCES(H, num_const, num_move, num_copy, num_inst);
        InstanceCounter<H>::update(num_const, num_move, num_copy, num_inst);

        EXPECT_EQ(3.0, h_fct4.call({ 3.0, 2.0 }));
    }
    EXPECT_EQ(num_const, InstanceCounter<H>::getNumberOfConstructions());
    EXPECT_EQ(num_copy, InstanceCounter<H>::getNumberOfCopies());
    // If zero instances are left, the destructor has been called the right
    // number of times, i.e., all internal casts in NamedFunction have been
    // successful.
    EXPECT_EQ(0, InstanceCounter<H>::getNumberOfInstances());
}
