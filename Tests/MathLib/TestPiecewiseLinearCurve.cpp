/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   TestPiecewiseLinearCurve.cpp
 *
 * Created on November 11, 2016, 10:58 AM
 */

#include <cmath>
#include <limits>

#include "gtest/gtest.h"
#include "TestTools.h"

#include "BaseLib/ConfigTree.h"
#include "MathLib/Curve/PiecewiseLinearCurve.h"
#include "MathLib/Curve/CreatePiecewiseLinearCurve.h"

std::unique_ptr<MathLib::PiecewiseLinearCurve> createPiecewiseLinearCurve(
    const char xml[])
{
    auto const ptree = readXml(xml);
    BaseLib::ConfigTree conf(ptree, "", BaseLib::ConfigTree::onerror,
                             BaseLib::ConfigTree::onwarning);
    auto const& sub_config = conf.getConfigSubtree("curve");
    return MathLib::createPiecewiseLinearCurve(sub_config);
}

TEST(MathLibCurve, PiecewiseLinearCurveParsing)
{
    const char xml[] =
        "<curve>"
        "   <type>PiecewiseLinear</type>"
        "   <coords> 0.2 0.4 0.5 0.6 0.7</coords>"
        "   <values> 20  10. 5.  3.  2. </values> "
        "</curve>";
    auto const curve = createPiecewiseLinearCurve(xml);

    std::vector<double> x = {0.2, 0.4, 0.5, 0.6, 0.7};
    std::vector<double> y = {
        20, 10, 5., 3., 2,
    };

    ASSERT_EQ(true, curve->isMonotonic());

    // Get inverse values and compare them
    for (std::size_t i = 0; i < x.size(); ++i)
    {
        ASSERT_NEAR(y[i], curve->getValue(x[i]),
                    std::numeric_limits<double>::epsilon());
        ASSERT_NEAR(x[i], curve->getVariable(y[i]),
                    std::numeric_limits<double>::epsilon());
    }
}

TEST(MathLibCurve, MonotonicIncreasePiecewiseLinearCurve)
{
    const std::size_t size = 100;
    std::vector<double> x;
    x.reserve(size);
    std::vector<double> y;
    y.reserve(size);

    const double dx = 10.0 / size;
    for (std::size_t i = 0; i < size; ++i)
    {
        const double xi = i * dx;
        x.push_back(xi);
        y.push_back(std::exp(xi));
    }

    std::vector<double> x_cpy = x;
    std::vector<double> y_cpy = y;

    MathLib::PiecewiseLinearCurve curve =
        MathLib::PiecewiseLinearCurve(std::move(x_cpy), std::move(y_cpy));

    ASSERT_EQ(true, curve.isMonotonic());

    // Get inverse values and compare them
    for (std::size_t i = 0; i < size; ++i)
    {
        ASSERT_NEAR(x[i], curve.getVariable(y[i]),
                    std::numeric_limits<double>::epsilon());
    }
}

TEST(MathLibCurve, MonotonicDecreasePiecewiseLinearCurve)
{
    const std::size_t size = 100;
    std::vector<double> x;
    x.reserve(size);
    std::vector<double> y;
    y.reserve(size);

    const double dx = 10.0 / size;
    for (std::size_t i = 0; i < size; ++i)
    {
        const double xi = i * dx;
        x.push_back(xi);
        y.push_back(std::exp(-xi));
    }

    std::vector<double> x_cpy = x;
    std::vector<double> y_cpy = y;

    MathLib::PiecewiseLinearCurve curve =
        MathLib::PiecewiseLinearCurve(std::move(x_cpy), std::move(y_cpy));

    ASSERT_EQ(true, curve.isMonotonic());

    // Get inverse values and compare them
    for (std::size_t i = 0; i < size; ++i)
    {
        ASSERT_NEAR(x[i], curve.getVariable(y[i]),
                    std::numeric_limits<double>::epsilon());
    }
}
