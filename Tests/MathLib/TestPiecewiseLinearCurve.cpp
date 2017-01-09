/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
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
#include "Tests/TestTools.h"

#include "BaseLib/ConfigTree.h"

#include "MathLib/Curve/CreatePiecewiseLinearCurve.h"
#include "MathLib/Curve/PiecewiseLinearMonotonicCurve.h"

template <typename CurveType>
std::unique_ptr<CurveType> createPiecewiseLinearCurve(const char xml[])
{
    auto const ptree = readXml(xml);
    BaseLib::ConfigTree conf(ptree, "", BaseLib::ConfigTree::onerror,
                             BaseLib::ConfigTree::onwarning);
    auto const& sub_config = conf.getConfigSubtree("curve");
    return MathLib::createPiecewiseLinearCurve<CurveType>(sub_config);
}

TEST(MathLibCurve, PiecewiseLinearCurveParsing)
{
    const char xml[] =
        "<curve>"
        "   <coords> 0.2 0.4 0.5 0.6 0.7</coords>"
        "   <values> 20  10. 5.  3.  2. </values> "
        "</curve>";
    auto const curve =
        createPiecewiseLinearCurve<MathLib::PiecewiseLinearMonotonicCurve>(xml);

    std::vector<double> x = {0.2, 0.4, 0.5, 0.6, 0.7};
    std::vector<double> y = {
        20, 10, 5., 3., 2,
    };

    // Get inverse values and compare them
    for (std::size_t i = 0; i < x.size(); ++i)
    {
        ASSERT_NEAR(y[i], curve->getValue(x[i]),
                    std::numeric_limits<double>::epsilon());
        ASSERT_NEAR(x[i], curve->getInversVariable(y[i]),
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

    MathLib::PiecewiseLinearMonotonicCurve curve =
        MathLib::PiecewiseLinearMonotonicCurve(std::move(x_cpy),
                                               std::move(y_cpy));

    // Get inverse values and compare them
    for (std::size_t i = 0; i < size; ++i)
    {
        ASSERT_NEAR(x[i], curve.getInversVariable(y[i]),
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

    MathLib::PiecewiseLinearMonotonicCurve curve =
        MathLib::PiecewiseLinearMonotonicCurve(std::move(x_cpy),
                                               std::move(y_cpy));

    // Get inverse values and compare them
    for (std::size_t i = 0; i < size; ++i)
    {
        ASSERT_NEAR(x[i], curve.getInversVariable(y[i]),
                    std::numeric_limits<double>::epsilon());
    }
}
