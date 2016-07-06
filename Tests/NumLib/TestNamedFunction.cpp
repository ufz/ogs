/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
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

    DBUG("calling %s", caller.getCallExpression("g").c_str());
    EXPECT_EQ(g_inst.g(x), caller.call("g", {x, y}));

    auto const f_caller = caller.getSpecialFunction("f");
    DBUG("calling %s", caller.getCallExpression("f").c_str());
    EXPECT_EQ(f_inst.f(g_inst.g(x), y), f_caller.call({x, y}));
}

TEST(NumLib, NamedFunctionNoLeaks)
{
    {
        H h_inst(1.0);

        auto h_fct = NumLib::NamedFunction("h", {"x", "y"},
                                           BaseLib::easyBind(&H::h, h_inst));
        EXPECT_EQ(5.0, h_fct.call({ 2.0, 3.0 }));
        h_inst.setZ(2.0);
        EXPECT_EQ(4.0, h_fct.call({ 2.0, 3.0 }));

        // Pass a temporary H object. NamedFunction will implicitly do memory
        // management.
        auto h_fct2 = NumLib::NamedFunction("h", {"x", "y"},
                                            BaseLib::easyBind(&H::h, H{3.0}));
        EXPECT_EQ(5.0, h_fct2.call({ 2.0, 4.0 }));

        // copy
        auto h_fct3 = h_fct2;
        EXPECT_EQ(-1.0, h_fct3.call({ 2.0, 1.0 }));

        // move
        auto h_fct4 = std::move(h_fct3);
        EXPECT_EQ(3.0, h_fct4.call({ 3.0, 2.0 }));
    }
    EXPECT_EQ(2, InstanceCounter<H>::getNumberOfConstructions());
    EXPECT_EQ(1, InstanceCounter<H>::getNumberOfCopies());
    // If zero instances are left, the destructor has been called the right
    // number of times, i.e., all internal casts in NamedFunction have been
    // successful.
    EXPECT_EQ(0, InstanceCounter<H>::getNumberOfInstances());
}
