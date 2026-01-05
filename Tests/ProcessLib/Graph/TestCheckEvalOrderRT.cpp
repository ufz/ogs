// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include "ProcessLib/Graph/CheckEvalOrderRT.h"

namespace TestCheckEvalOrderRT
{
struct NoOp
{
    void eval();
};

struct OutD
{
    void eval(double&) {}
};

struct InD
{
    void eval(double const&) {}
};

struct InSOutD
{
    void eval(double&, std::string const&) {}
};

struct OutCD
{
    void eval(double&, char&) {}
};

struct InCSOutIU
{
    void eval(char const&, unsigned&, std::string const&, int&);
};

struct InIU
{
    void eval(unsigned const&, int const&);
};
}  // namespace TestCheckEvalOrderRT

TEST(ProcessLib_Graph_CheckEvalOrderRT, Success)
{
    using namespace TestCheckEvalOrderRT;

    {
        using Models = std::tuple<NoOp>;
        using Inputs = std::tuple<>;

        ASSERT_TRUE(
            (ProcessLib::Graph::isEvalOrderCorrectRT<Models, Inputs>()));
    }

    {
        using Models = std::tuple<OutD>;
        using Inputs = std::tuple<>;

        ASSERT_TRUE(
            (ProcessLib::Graph::isEvalOrderCorrectRT<Models, Inputs>()));
    }

    {
        using Models = std::tuple<InD>;
        using Inputs = std::tuple<double>;

        ASSERT_TRUE(
            (ProcessLib::Graph::isEvalOrderCorrectRT<Models, Inputs>()));
    }

    {
        using Models = std::tuple<OutD, InD>;
        using Inputs = std::tuple<>;

        ASSERT_TRUE(
            (ProcessLib::Graph::isEvalOrderCorrectRT<Models, Inputs>()));
    }

    {
        using Models = std::tuple<InSOutD, NoOp, InD>;
        using Inputs = std::tuple<std::string>;

        ASSERT_TRUE(
            (ProcessLib::Graph::isEvalOrderCorrectRT<Models, Inputs>()));
    }

    {
        using Models = std::tuple<OutCD, NoOp, InCSOutIU, InD, InIU>;
        using Inputs = std::tuple<std::string>;

        ASSERT_TRUE(
            (ProcessLib::Graph::isEvalOrderCorrectRT<Models, Inputs>()));
    }

    {
        using Models = std::tuple<NoOp, OutCD, InD, InCSOutIU, InIU>;
        using Inputs = std::tuple<std::string>;

        ASSERT_TRUE(
            (ProcessLib::Graph::isEvalOrderCorrectRT<Models, Inputs>()));
    }
}

TEST(ProcessLib_Graph_CheckEvalOrderRT, Fail)
{
    using namespace TestCheckEvalOrderRT;

    {
        using Models = std::tuple<InD>;
        using Inputs = std::tuple<>;

        ASSERT_FALSE(
            (ProcessLib::Graph::isEvalOrderCorrectRT<Models, Inputs>()))
            << "Expected to fail since the required double value is not "
               "provided as an input.";
    }

    {
        using Models = std::tuple<InD, OutD>;
        using Inputs = std::tuple<>;

        ASSERT_FALSE(
            (ProcessLib::Graph::isEvalOrderCorrectRT<Models, Inputs>()))
            << "Expected to fail: double value is consumed before being "
               "computed.";
    }

    {
        using Models = std::tuple<NoOp, InSOutD, OutD, InD>;
        using Inputs = std::tuple<std::string>;

        ASSERT_FALSE(
            (ProcessLib::Graph::isEvalOrderCorrectRT<Models, Inputs>()))
            << "Expected to fail: double value has two producers.";
    }

    {
        using Models = std::tuple<OutCD, InIU, NoOp, InCSOutIU, InD>;
        using Inputs = std::tuple<std::string>;

        ASSERT_FALSE(
            (ProcessLib::Graph::isEvalOrderCorrectRT<Models, Inputs>()))
            << "Expected to fail: wrong evaluation order.";
    }

    {
        using Models = std::tuple<NoOp, OutCD, InD, OutD, InCSOutIU, InIU>;
        using Inputs = std::tuple<std::string>;

        ASSERT_FALSE(
            (ProcessLib::Graph::isEvalOrderCorrectRT<Models, Inputs>()))
            << "Expected to fail: double value has two producers.";
    }
}
