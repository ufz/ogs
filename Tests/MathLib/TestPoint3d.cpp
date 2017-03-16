/**
 * @file TestPoint3d.cpp
 * @author Thomas Fischer
 * @date Nov 8, 2012
 *
 * @copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <ctime>
#include <gtest/gtest.h>
#include <autocheck/autocheck.hpp>

#include "MathLib/Point3d.h"

#include "Tests/AutoCheckTools.h"

using namespace MathLib;
namespace ac = autocheck;

struct MathLibPoint3d : public ::testing::Test
{
    ac::randomTupleGenerator<double, 3> tupleGen;
    ac::cons_generator<MathLib::Point3d, ac::randomTupleGenerator<double, 3>>
        pointGenerator{tupleGen};

    ac::randomCoordinateIndexGenerator<unsigned, 3>
        coordGenerator;  // any of {0, 1, 2}
    ac::gtest_reporter gtest_reporter;
};

TEST_F(MathLibPoint3d, ComparisonOperatorLessEqSamePoint)
{
    // A point is always less or equal to itself and its copy.
    auto samePointLessEqualCompare = [](MathLib::Point3d const& p)
    {
        const auto& q = p;
        return lessEq(p, p) && lessEq(p, q) && lessEq(q, p);
    };

    ac::check<MathLib::Point3d>(samePointLessEqualCompare, 1000,
                                ac::make_arbitrary(pointGenerator),
                                gtest_reporter);
}

TEST_F(MathLibPoint3d, ComparisonOperatorLessEqualLargePerturbation)
{
    // A point with any big, positive value added to one of its coordinates is
    // never smaller or equal to the original point.
    // And the original point is always smaller or equal to the perturbed point.
    auto pointWithLargeAddedValue =
        [](MathLib::Point3d const& p, double const perturbation,
           unsigned const coordinate)
    {
        auto q = p;
        q[coordinate] = q[coordinate] + perturbation;
        return !lessEq(q, p) && lessEq(p, q);
    };

    auto eps = std::numeric_limits<double>::epsilon();

    ac::check<MathLib::Point3d, double, unsigned>(
        pointWithLargeAddedValue, 10000,
        ac::make_arbitrary(pointGenerator,
                           ac::map(&ac::absoluteValue, ac::generator<double>()),
                           coordGenerator)
            .discard_if(
                [&eps](MathLib::Point3d const&, double const v, unsigned const)
                {
                    return !(v > eps);
                }),
        gtest_reporter);
}

TEST_F(MathLibPoint3d, ComparisonOperatorLessEqualSmallPerturbation)
{
    // A point with any non-zero value smaller than epsilon/2 added to one of
    // its
    // coordinates is always less or equal to the original point.
    auto pointWithSmallAddedValue =
        [](MathLib::Point3d const& p, double const perturbation,
           unsigned const coordinate)
    {
        auto q = p;
        q[coordinate] = q[coordinate] + perturbation;
        return lessEq(p, q) && lessEq(q, p);
    };

    auto eps = std::numeric_limits<double>::epsilon();

    ac::check<MathLib::Point3d, double, unsigned>(
        pointWithSmallAddedValue, 10000,
        ac::make_arbitrary(pointGenerator,
                           ac::progressivelySmallerGenerator<double>(eps / 2),
                           coordGenerator),
        gtest_reporter);
}

TEST_F(MathLibPoint3d, ComparisonOperatorEqualSamePoint)
{
    // A point is always equal to itself and its copy.
    auto samePointEqualCompare = [](MathLib::Point3d const& p)
    {
        const auto& q = p;
        return (p == p) && (p == q) && (q == p);
    };

    ac::check<MathLib::Point3d>(samePointEqualCompare, 100,
                                ac::make_arbitrary(pointGenerator),
                                gtest_reporter);
}

TEST_F(MathLibPoint3d, ComparisonOperatorEqualLargePerturbation)
{
    // A point with any big, non-zero value added to one of its coordinates is
    // never equal to the original point.
    auto pointWithLargeAddedValue =
        [](MathLib::Point3d const& p, double const perturbation,
           unsigned const coordinate)
    {
        auto q = p;
        q[coordinate] = q[coordinate] + perturbation;
        return !(p == q) && !(q == p);
    };

    auto eps = std::numeric_limits<double>::epsilon();

    ac::check<MathLib::Point3d, double, unsigned>(
        pointWithLargeAddedValue, 10000,
        ac::make_arbitrary(pointGenerator, ac::generator<double>(),
                           coordGenerator)
            .discard_if(
                [&eps](MathLib::Point3d const&, double const v, unsigned const)
                {
                    return !(v > eps);
                }),
        gtest_reporter);
}

TEST_F(MathLibPoint3d, ComparisonOperatorEqualSmallPerturbation)
{
    // A point with any non-zero value smaller than epsilon/2 added to one of
    // its
    // coordinates is always equal to the original point.
    auto pointWithSmallAddedValue =
        [](MathLib::Point3d const& p, double const perturbation,
           unsigned const coordinate)
    {
        auto q = p;
        q[coordinate] = q[coordinate] + perturbation;
        return (p == q) && (q == p);
    };

    auto eps = std::numeric_limits<double>::epsilon();

    ac::check<MathLib::Point3d, double, unsigned>(
        pointWithSmallAddedValue, 1000,
        ac::make_arbitrary(pointGenerator,
                           ac::progressivelySmallerGenerator<double>(eps / 2),
                           coordGenerator),
        gtest_reporter);
}

// test for operator<
TEST_F(MathLibPoint3d, ComparisonOperatorLessSamePoint)
{
    // A point is never less than itself or its copy.
    auto samePointLessCompare = [](MathLib::Point3d const& p)
    {
        const auto& q = p;
        return !(p < p) && !(p < q) && !(q < p);
    };

    ac::check<MathLib::Point3d>(samePointLessCompare, 100,
                                ac::make_arbitrary(pointGenerator),
                                gtest_reporter);
}

TEST_F(MathLibPoint3d, ComparisonOperatorLessLargePerturbation)
{
    // A point with any positive value added to one of its coordinates is
    // always larger then the original point.
    auto pointWithAddedValue = [](MathLib::Point3d const& p, double const eps,
                                  unsigned const coordinate)
    {
        auto q = p;
        q[coordinate] = q[coordinate] + eps;
        return (p < q) && !(q < p);
    };

    ac::check<MathLib::Point3d, double, unsigned>(
        pointWithAddedValue, 1000,
        ac::make_arbitrary(pointGenerator,
                           ac::map(&ac::absoluteValue, ac::generator<double>()),
                           coordGenerator)
            .discard_if(
                [](MathLib::Point3d const&, double const eps, unsigned const)
                {
                    return eps == 0;
                }),
        gtest_reporter);
}

TEST_F(MathLibPoint3d, ComparisonOperatorLessSmallPerturbation)
{
    // A point with any positive value subtracted from one of its coordinates is
    // always smaller then the original point.
    auto pointWithSubtractedValue = [](
        MathLib::Point3d const& p, double const eps, unsigned const coordinate)
    {
        auto q = p;
        q[coordinate] = q[coordinate] - eps;
        return (q < p) && !(p < q);
    };

    ac::check<MathLib::Point3d, double, unsigned>(
        pointWithSubtractedValue, 1000,
        ac::make_arbitrary(pointGenerator,
                           ac::map(&ac::absoluteValue, ac::generator<double>()),
                           coordGenerator)
            .discard_if(
                [](MathLib::Point3d const&, double const eps, unsigned const)
                {
                    return eps == 0;
                }),
        gtest_reporter);
}
