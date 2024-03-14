/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <gtest/gtest.h>

#include <ctime>
#include <memory>
#include <random>

#include "BaseLib/Algorithm.h"
#include "GeoLib/AABB.h"
#include "GeoLib/OctTree.h"
#include "GeoLib/Point.h"

namespace
{
Eigen::Vector3d convertToEigen(MathLib::Point3d point3d)
{
    return Eigen::Vector3d{point3d[0], point3d[1], point3d[2]};
}
}  // namespace

class GeoLibOctTree : public testing::Test
{
public:
    using VectorOfPoints = std::vector<GeoLib::Point*>;

    GeoLibOctTree() = default;
    ~GeoLibOctTree() override { BaseLib::cleanupVectorElements(ps_ptr); }

#ifndef NDEBUG
    template <std::size_t MAX_POINTS>
    void checkOctTreeChildsNonNullptr(
        GeoLib::OctTree<GeoLib::Point, MAX_POINTS> const& oct_tree) const
    {
        ASSERT_NE(nullptr, oct_tree.getChild(0));
        ASSERT_NE(nullptr, oct_tree.getChild(1));
        ASSERT_NE(nullptr, oct_tree.getChild(2));
        ASSERT_NE(nullptr, oct_tree.getChild(3));
        ASSERT_NE(nullptr, oct_tree.getChild(4));
        ASSERT_NE(nullptr, oct_tree.getChild(5));
        ASSERT_NE(nullptr, oct_tree.getChild(6));
        ASSERT_NE(nullptr, oct_tree.getChild(7));
    }
#endif

#ifndef NDEBUG
    template <std::size_t MAX_POINTS>
    void checkOctTreeChildsNullptr(
        GeoLib::OctTree<GeoLib::Point, MAX_POINTS> const& oct_tree) const
    {
        ASSERT_EQ(nullptr, oct_tree.getChild(0));
        ASSERT_EQ(nullptr, oct_tree.getChild(1));
        ASSERT_EQ(nullptr, oct_tree.getChild(2));
        ASSERT_EQ(nullptr, oct_tree.getChild(3));
        ASSERT_EQ(nullptr, oct_tree.getChild(4));
        ASSERT_EQ(nullptr, oct_tree.getChild(5));
        ASSERT_EQ(nullptr, oct_tree.getChild(6));
        ASSERT_EQ(nullptr, oct_tree.getChild(7));
    }
#endif

protected:
    void generateEquidistantPoints3d(std::size_t const n = 11)
    {
        for (std::size_t k(0); k < n; ++k)
        {
            double const z(k - (n - 1) / 2.0);
            for (std::size_t j(0); j < n; ++j)
            {
                double const y(j - (n - 1) / 2.0);
                for (std::size_t i(0); i < n; ++i)
                {
                    ps_ptr.push_back(
                        new GeoLib::Point(i - (n - 1) / 2.0, y, z));
                }
            }
        }
    }

    void generateEquidistantPoints3dUnitCube(std::size_t const n = 11)
    {
        for (std::size_t k(0); k < n; ++k)
        {
            double const z(k / (n - 1.0));
            for (std::size_t j(0); j < n; ++j)
            {
                double const y(j / (n - 1.0));
                for (std::size_t i(0); i < n; ++i)
                {
                    ps_ptr.push_back(new GeoLib::Point(i / (n - 1.0), y, z));
                }
            }
        }
    }

    void generateEquidistantPointsUnitSquare(std::size_t const n)
    {
        for (std::size_t i = 0; i <= n; ++i)
        {
            for (std::size_t j = 0; j <= n; ++j)
            {
                ps_ptr.push_back(new GeoLib::Point(i / double(n), j / double(n),
                                                   0, i * (n + 1) + j));
            }
        }
    }

    std::tuple<Eigen::Vector3d const, Eigen::Vector3d const,
               std::unique_ptr<GeoLib::OctTree<GeoLib::Point, 2>>>
    generateOctTreeFromPointSet(double const eps)
    {
        GeoLib::AABB aabb(ps_ptr.begin(), ps_ptr.end());
        auto const& min(aabb.getMinPoint());
        auto const& max(aabb.getMaxPoint());
        std::unique_ptr<GeoLib::OctTree<GeoLib::Point, 2>> oct_tree(
            GeoLib::OctTree<GeoLib::Point, 2>::createOctTree(min, max, eps));
        for (auto* p : ps_ptr)
        {
            GeoLib::Point* ret_pnt(nullptr);
            // the insertion of points into the OctTree is already tested
            if (!oct_tree->addPoint(p, ret_pnt))
            {
                OGS_FATAL("Could not insert point into OctTree");
            }
        }
        return std::tuple<Eigen::Vector3d const, Eigen::Vector3d const,
                          std::unique_ptr<GeoLib::OctTree<GeoLib::Point, 2>>>(
            min, max, std::move(oct_tree));
    }

protected:
    VectorOfPoints ps_ptr;
};

TEST_F(GeoLibOctTree, TestWithEquidistantPoints3d)
{
    generateEquidistantPoints3d();
    double const eps(10 * std::numeric_limits<double>::epsilon());
    std::unique_ptr<GeoLib::OctTree<GeoLib::Point, 2>> oct_tree(
        GeoLib::OctTree<GeoLib::Point, 2>::createOctTree(
            convertToEigen(*ps_ptr.front()), convertToEigen(*ps_ptr.back()),
            eps));

#ifndef NDEBUG
    Eigen::Vector3d const& ll(oct_tree->getLowerLeftCornerPoint());
    Eigen::Vector3d const& ur(oct_tree->getUpperRightCornerPoint());

    EXPECT_EQ((*ps_ptr.front())[0], ll[0]);
    EXPECT_EQ((*ps_ptr.front())[1], ll[1]);
    EXPECT_EQ((*ps_ptr.front())[2], ll[2]);

    EXPECT_NEAR((*ps_ptr.back())[0], ur[0], (ur[0] - ll[0]) * 1e-6);
    EXPECT_NEAR((*ps_ptr.back())[1], ur[1], (ur[1] - ll[1]) * 1e-6);
    EXPECT_NEAR((*ps_ptr.back())[2], ur[2], (ur[2] - ll[2]) * 1e-6);

    checkOctTreeChildsNullptr<2>(*oct_tree);

    ASSERT_EQ(static_cast<std::size_t>(0), oct_tree->getPointVector().size());
#endif

    GeoLib::Point* ret_pnt(nullptr);
    // insert the first point
    ASSERT_TRUE(oct_tree->addPoint(ps_ptr[0], ret_pnt));

    // make a range query
    MathLib::Point3d const min(
        std::array<double, 3>{{(*(ps_ptr[0]))[0] - eps, (*(ps_ptr[0]))[1] - eps,
                               (*(ps_ptr[0]))[2] - eps}});
    MathLib::Point3d const max(
        std::array<double, 3>{{(*(ps_ptr[0]))[0] + eps, (*(ps_ptr[0]))[1] + eps,
                               (*(ps_ptr[0]))[2] + eps}});
    std::vector<GeoLib::Point*> query_pnts;
    oct_tree->getPointsInRange(min, max, query_pnts);
    ASSERT_EQ(1u, query_pnts.size());

#ifndef NDEBUG
    ASSERT_EQ(static_cast<std::size_t>(1), oct_tree->getPointVector().size());
#endif

    // try to insert the first point a second time
    ASSERT_FALSE(oct_tree->addPoint(ps_ptr[0], ret_pnt));
#ifndef NDEBUG
    ASSERT_EQ(static_cast<std::size_t>(1), oct_tree->getPointVector().size());
#endif

    // insert the second point
    ASSERT_TRUE(oct_tree->addPoint(ps_ptr[1], ret_pnt));
#ifndef NDEBUG
    ASSERT_EQ(static_cast<std::size_t>(2), oct_tree->getPointVector().size());
    checkOctTreeChildsNullptr<2>(*oct_tree);
#endif

    // insert a third point -> there should be subtrees
    ASSERT_TRUE(oct_tree->addPoint(ps_ptr[2], ret_pnt));
#ifndef NDEBUG
    ASSERT_EQ(static_cast<std::size_t>(0), oct_tree->getPointVector().size());
    checkOctTreeChildsNonNullptr<2>(*oct_tree);

    // all inserted points are in the SWL -> there should be another subtree
    // level
    checkOctTreeChildsNonNullptr<2>(*(oct_tree->getChild(2)));

    // still all inserted points are in the SWL of the SWL
    // -> there should be another subtree level
    checkOctTreeChildsNonNullptr<2>(*(oct_tree->getChild(2)->getChild(2)));

    checkOctTreeChildsNullptr<2>(*(oct_tree->getChild(2)->getChild(0)));
    checkOctTreeChildsNullptr<2>(*(oct_tree->getChild(2)->getChild(1)));
    checkOctTreeChildsNullptr<2>(*(oct_tree->getChild(2)->getChild(3)));
    checkOctTreeChildsNullptr<2>(*(oct_tree->getChild(2)->getChild(4)));
    checkOctTreeChildsNullptr<2>(*(oct_tree->getChild(2)->getChild(5)));
    checkOctTreeChildsNullptr<2>(*(oct_tree->getChild(2)->getChild(6)));
    checkOctTreeChildsNullptr<2>(*(oct_tree->getChild(2)->getChild(7)));

    ASSERT_EQ(static_cast<std::size_t>(2),
              oct_tree->getChild(2)
                  ->getChild(2)
                  ->getChild(2)
                  ->getPointVector()
                  .size());
    ASSERT_EQ(static_cast<std::size_t>(1),
              oct_tree->getChild(2)
                  ->getChild(2)
                  ->getChild(3)
                  ->getPointVector()
                  .size());
#endif

    ASSERT_TRUE(oct_tree->addPoint(ps_ptr[3], ret_pnt));
#ifndef NDEBUG
    ASSERT_EQ(static_cast<std::size_t>(1),
              oct_tree->getChild(2)->getChild(3)->getPointVector().size());
#endif

    GeoLib::Point range_query_ll(*(ps_ptr.front()));
    GeoLib::Point range_query_ur(*(ps_ptr[ps_ptr.size() / 2]));
    std::vector<GeoLib::Point*> result;
    oct_tree->getPointsInRange(range_query_ll, range_query_ur, result);
    ASSERT_EQ(static_cast<std::size_t>(4), result.size());

    result.clear();
    range_query_ur[0] = -2.5;
    range_query_ur[1] = -2.5;
    range_query_ur[2] = -2.5;
    oct_tree->getPointsInRange(range_query_ll, range_query_ur, result);
    ASSERT_EQ(static_cast<std::size_t>(3), result.size());

    // insert some points not resulting in a further refinement of SWL
    for (std::size_t k(4); k < 11; ++k)
    {
        ASSERT_TRUE(oct_tree->addPoint(ps_ptr[k], ret_pnt));
    }

    result.clear();
    oct_tree->getPointsInRange(range_query_ll, range_query_ur, result);
    ASSERT_EQ(static_cast<std::size_t>(3), result.size());

    // insert some points *resulting* in a further refinement of SWL
    for (std::size_t k(11); k < 25; ++k)
    {
        ASSERT_TRUE(oct_tree->addPoint(ps_ptr[k], ret_pnt));
    }

    result.clear();
    oct_tree->getPointsInRange(range_query_ll, range_query_ur, result);
    ASSERT_EQ(static_cast<std::size_t>(9), result.size());

    // insert all points with z = -5.0 - this does not result in a further
    // refinement of SWL
    for (std::size_t k(25); k < 121; ++k)
    {
        ASSERT_TRUE(oct_tree->addPoint(ps_ptr[k], ret_pnt));
    }

    result.clear();
    oct_tree->getPointsInRange(range_query_ll, range_query_ur, result);
    ASSERT_EQ(static_cast<std::size_t>(9), result.size());

    result.clear();
    range_query_ur[0] = -3.75;
    range_query_ur[1] = -3.75;
    range_query_ur[2] = -3.75;
    oct_tree->getPointsInRange(range_query_ll, range_query_ur, result);
    ASSERT_EQ(static_cast<std::size_t>(4), result.size());

    result.clear();
    range_query_ur[0] = -4.25;
    range_query_ur[1] = -4.25;
    range_query_ur[2] = -4.25;
    oct_tree->getPointsInRange(range_query_ll, range_query_ur, result);
    ASSERT_EQ(static_cast<std::size_t>(1), result.size());

    result.clear();
    range_query_ll[0] = -4.75;
    range_query_ll[1] = -4.75;
    range_query_ll[2] = -4.75;
    oct_tree->getPointsInRange(range_query_ll, range_query_ur, result);
    for (auto p : result)
    {
        std::cout << *p << "\n";
    }
    ASSERT_EQ(static_cast<std::size_t>(0), result.size());

    result.clear();
    range_query_ll[0] = -5;
    range_query_ll[1] = -5;
    range_query_ll[2] = -5;
    range_query_ur[0] = -0.25;
    range_query_ur[1] = -4.75;
    range_query_ur[2] = -4.75;
    oct_tree->getPointsInRange(range_query_ll, range_query_ur, result);
    ASSERT_EQ(static_cast<std::size_t>(5), result.size());
}

TEST_F(GeoLibOctTree, TestWithAlternatingPoints3d)
{
    // this case is not correctly handled by lexicographical sorting
    double const eps(1e-1);
    double const small_displacement(1e-2);
    ps_ptr.push_back(new GeoLib::Point(0, 0, 0, 0));
    ps_ptr.push_back(new GeoLib::Point(2 * small_displacement, 0, 0, 1));
    ps_ptr.push_back(new GeoLib::Point(small_displacement, 1, 0, 2));
    ps_ptr.push_back(new GeoLib::Point(4 * small_displacement, 0, 0, 3));
    ps_ptr.push_back(new GeoLib::Point(3 * small_displacement, 1, 0, 4));
    ps_ptr.push_back(new GeoLib::Point(6 * small_displacement, 0, 0, 5));
    ps_ptr.push_back(new GeoLib::Point(5 * small_displacement, 1, 0, 6));

    GeoLib::AABB const aabb(ps_ptr.cbegin(), ps_ptr.cend());
    auto const& min(aabb.getMinPoint());
    auto const& max(aabb.getMaxPoint());
    std::unique_ptr<GeoLib::OctTree<GeoLib::Point, 8>> oct_tree(
        GeoLib::OctTree<GeoLib::Point, 8>::createOctTree(min, max, eps));

    // pt_ptr[0] should be inserted correctly
    GeoLib::Point* ret_pnt(nullptr);
    ASSERT_TRUE(oct_tree->addPoint(ps_ptr[0], ret_pnt));
    ASSERT_EQ(ps_ptr[0], ret_pnt);
    // ps_ptr[1] is in the eps-environment of ps_ptr[0]
    ret_pnt = nullptr;
    ASSERT_FALSE(oct_tree->addPoint(ps_ptr[1], ret_pnt));
    ASSERT_EQ(ps_ptr[0], ret_pnt);
    // pt_ptr[2] should be inserted correctly
    ret_pnt = nullptr;
    ASSERT_TRUE(oct_tree->addPoint(ps_ptr[2], ret_pnt));
    ASSERT_EQ(ps_ptr[2], ret_pnt);
    // ps_ptr[3] is in the eps-environment of ps_ptr[0]
    ret_pnt = nullptr;
    ASSERT_FALSE(oct_tree->addPoint(ps_ptr[3], ret_pnt));
    ASSERT_EQ(ps_ptr[0], ret_pnt);
    // ps_ptr[4] is in the eps-environment of ps_ptr[2]
    ret_pnt = nullptr;
    ASSERT_FALSE(oct_tree->addPoint(ps_ptr[4], ret_pnt));
    ASSERT_EQ(ps_ptr[2], ret_pnt);
    // ps_ptr[5] is in the eps-environment of ps_ptr[0]
    ret_pnt = nullptr;
    ASSERT_FALSE(oct_tree->addPoint(ps_ptr[5], ret_pnt));
    ASSERT_EQ(ps_ptr[0], ret_pnt);
    // ps_ptr[6] is in the eps-environment of ps_ptr[2]
    ret_pnt = nullptr;
    ASSERT_FALSE(oct_tree->addPoint(ps_ptr[6], ret_pnt));
    ASSERT_EQ(ps_ptr[2], ret_pnt);
}

TEST_F(GeoLibOctTree, TestRangeQueryOnCube)
{
    generateEquidistantPoints3d(21);

    // create OctTree
    GeoLib::AABB const aabb(ps_ptr.cbegin(), ps_ptr.cend());
    auto const& aabb_min(aabb.getMinPoint());
    auto const& aabb_max(aabb.getMaxPoint());
    double const eps(0.5);
    std::unique_ptr<GeoLib::OctTree<GeoLib::Point, 2>> oct_tree(
        GeoLib::OctTree<GeoLib::Point, 2>::createOctTree(aabb_min, aabb_max,
                                                         eps));
    // fill OctTree with the test points
    for (auto p : ps_ptr)
    {
        GeoLib::Point* ret_pnt(nullptr);
        ASSERT_TRUE(oct_tree->addPoint(p, ret_pnt));
        ASSERT_EQ(p, ret_pnt);
    }

    // do a range query for every test point - only the test point should be
    // found
    for (auto const* point : ps_ptr)
    {
        std::vector<GeoLib::Point*> found_points;
        Eigen::Vector3d min = point->asEigenVector3d().array() - eps;
        Eigen::Vector3d max = point->asEigenVector3d().array() + eps;
        oct_tree->getPointsInRange(min, max, found_points);
        ASSERT_EQ(1u, found_points.size());
        ASSERT_EQ((*point)[0], (*found_points[0])[0]);
        ASSERT_EQ((*point)[1], (*found_points[0])[1]);
        ASSERT_EQ((*point)[2], (*found_points[0])[2]);
        ASSERT_EQ(point->getID(), found_points[0]->getID());
    }
}

TEST_F(GeoLibOctTree, TestOctTreeWithTwoEqualPoints)
{
    ps_ptr.push_back(new GeoLib::Point(0, 0, 0, 0));
    ps_ptr.push_back(new GeoLib::Point(0, 0, 0, 1));
    double const eps(0.0);

    GeoLib::AABB aabb(ps_ptr.begin(), ps_ptr.end());
    auto const& min(aabb.getMinPoint());
    auto const& max(aabb.getMaxPoint());
    std::unique_ptr<GeoLib::OctTree<GeoLib::Point, 2>> oct_tree(
        GeoLib::OctTree<GeoLib::Point, 2>::createOctTree(min, max, eps));

    GeoLib::Point* ret_pnt(nullptr);
    ASSERT_TRUE(oct_tree->addPoint(ps_ptr[0], ret_pnt));
    ASSERT_EQ(ps_ptr[0], ret_pnt);
    ASSERT_FALSE(oct_tree->addPoint(ps_ptr[1], ret_pnt));
    ASSERT_EQ(ps_ptr[0], ret_pnt);
}

TEST_F(GeoLibOctTree, TestOctTreeWithTwoEqualPointsOne)
{
    ps_ptr.push_back(new GeoLib::Point(1, 1, 1, 0));
    ps_ptr.push_back(new GeoLib::Point(1, 1, 1, 1));
    double const eps(0.0);

    GeoLib::AABB aabb(ps_ptr.begin(), ps_ptr.end());
    auto const& min(aabb.getMinPoint());
    auto const& max(aabb.getMaxPoint());
    std::unique_ptr<GeoLib::OctTree<GeoLib::Point, 2>> oct_tree(
        GeoLib::OctTree<GeoLib::Point, 2>::createOctTree(min, max, eps));

    GeoLib::Point* ret_pnt(nullptr);
    ASSERT_TRUE(oct_tree->addPoint(ps_ptr[0], ret_pnt));
    ASSERT_EQ(ps_ptr[0], ret_pnt);
    ASSERT_FALSE(oct_tree->addPoint(ps_ptr[1], ret_pnt));
    ASSERT_EQ(ps_ptr[0], ret_pnt);
}

TEST_F(GeoLibOctTree, TestOctTreeOnCubicDomain)
{
    ps_ptr.push_back(new GeoLib::Point(-1, -1, -1, 0));
    ps_ptr.push_back(new GeoLib::Point(1, 1, 1, 1));
    double const eps(0.0);

    GeoLib::AABB aabb(ps_ptr.begin(), ps_ptr.end());
    auto const& min(aabb.getMinPoint());
    auto const& max(aabb.getMaxPoint());
    std::unique_ptr<GeoLib::OctTree<GeoLib::Point, 2>> oct_tree(
        GeoLib::OctTree<GeoLib::Point, 2>::createOctTree(min, max, eps));

    GeoLib::Point* ret_pnt(nullptr);
    ASSERT_TRUE(oct_tree->addPoint(ps_ptr[0], ret_pnt));
    ASSERT_EQ(ps_ptr[0], ret_pnt);
    ASSERT_TRUE(oct_tree->addPoint(ps_ptr[1], ret_pnt));
    ASSERT_EQ(ps_ptr[1], ret_pnt);
}

TEST_F(GeoLibOctTree, TestAddPointOnSquareDomain)
{
    generateEquidistantPointsUnitSquare(3);
    double const eps = std::numeric_limits<double>::epsilon() * 0.5;
    GeoLib::AABB aabb(ps_ptr.begin(), ps_ptr.end());
    auto const& min(aabb.getMinPoint());
    auto const& max(aabb.getMaxPoint());
    std::unique_ptr<GeoLib::OctTree<GeoLib::Point, 2>> oct_tree(
        GeoLib::OctTree<GeoLib::Point, 2>::createOctTree(min, max, eps));
    for (auto* p : ps_ptr)
    {
        GeoLib::Point* ret_pnt(nullptr);
        ASSERT_TRUE(oct_tree->addPoint(p, ret_pnt));
        ASSERT_EQ(p, ret_pnt);
    }

    VectorOfPoints points_with_same_coordinates;
    for (auto const* point : ps_ptr)
    {
        points_with_same_coordinates.push_back(new GeoLib::Point(*point));
    }

    for (auto* p : points_with_same_coordinates)
    {
        GeoLib::Point* ret_pnt(nullptr);
        ASSERT_FALSE(oct_tree->addPoint(p, ret_pnt));
        ASSERT_EQ((*p)[0], (*ret_pnt)[0]);
        ASSERT_EQ((*p)[1], (*ret_pnt)[1]);
        ASSERT_EQ((*p)[2], (*ret_pnt)[2]);
        ASSERT_FALSE(p == ret_pnt);
    }
    BaseLib::cleanupVectorElements(points_with_same_coordinates);
}

TEST_F(GeoLibOctTree, TestRangeQueryForEntireDomain)
{
    std::size_t const n = 3;
    generateEquidistantPointsUnitSquare(n);
    double const eps = std::numeric_limits<double>::epsilon() * 0.5;
    auto [min, max, oct_tree] = generateOctTreeFromPointSet(eps);

    std::vector<GeoLib::Point*> query_points;
    // min and max from aabb -> all inserted points should be in query_pnts
    oct_tree->getPointsInRange(min, max, query_points);
    ASSERT_EQ((n + 1) * (n + 1), query_points.size());
}

TEST_F(GeoLibOctTree, TestRangeQueryEmptyRange)
{
    generateEquidistantPoints3dUnitCube(21);
    double const eps = std::numeric_limits<double>::epsilon() * 0.5;
    auto [min, max, oct_tree] = generateOctTreeFromPointSet(eps);

    for (auto const* point : ps_ptr)
    {
        std::vector<GeoLib::Point*> query_points;
        Eigen::Vector3d const min_p(point->asEigenVector3d());
        Eigen::Vector3d const max_p = min_p;
        oct_tree->getPointsInRange(min_p, max_p, query_points);
        ASSERT_EQ(0u, query_points.size());
    }
}

TEST_F(GeoLibOctTree, TestRangeQueryWithOutsideRange)
{
    generateEquidistantPoints3dUnitCube(21);
    double const eps = std::numeric_limits<double>::epsilon() * 0.5;
    auto [min, max, oct_tree] = generateOctTreeFromPointSet(eps);

    // range query for range outside the cube domain [min, max)
    std::vector<GeoLib::Point*> query_points;
    Eigen::Vector3d const min_p(max);
    Eigen::Vector3d const max_p = min_p.array() + 1.0;
    oct_tree->getPointsInRange(min_p, max_p, query_points);
    ASSERT_EQ(0u, query_points.size());
}
