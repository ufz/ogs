/**
 * \file
 * \author Thomas Fischer
 * \date   Oct 24, 2012
 * \brief  Tests of AABB class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <cstdlib>
#include <ctime>
#include <list>
#include <random>

#include "gtest/gtest.h"

#include "GeoLib/AABB.h"
#include "GeoLib/Point.h"
#include "MathLib/Point3d.h"

TEST(GeoLibAABB, RandomNumberOfPointersToRandomPoints)
{
    /* initialize random seed: */
     srand ( static_cast<unsigned>(time(nullptr)) );
     int n (rand() % 100000);
     int box_size (rand());
     double half_box_size(box_size/2);
     double minus_half_box_size(-half_box_size);

     // fill list with points
     std::list<GeoLib::Point*> pnts_list;
     for (int k(0); k<n; k++) {
         pnts_list.push_back(new GeoLib::Point(rand() % box_size - half_box_size, rand() % box_size - half_box_size, rand() % box_size - half_box_size));
     }

     // construct from list points a axis algined bounding box
     GeoLib::AABB aabb(pnts_list.begin(), pnts_list.end());

     MathLib::Point3d const& min_pnt(aabb.getMinPoint());
     MathLib::Point3d const& max_pnt(aabb.getMaxPoint());

     ASSERT_LE(minus_half_box_size, min_pnt[0]) << "coordinate 0 of min_pnt is smaller than " << minus_half_box_size;
     ASSERT_LE(minus_half_box_size, min_pnt[1]) << "coordinate 1 of min_pnt is smaller than " << minus_half_box_size;
     ASSERT_LE(minus_half_box_size, min_pnt[2]) << "coordinate 2 of min_pnt is smaller than " << minus_half_box_size;

    // since the interval is half-open we have to move the upper bound also
     half_box_size = std::nexttoward(half_box_size, std::numeric_limits<double>::max());
     ASSERT_GE(half_box_size, max_pnt[0]) << "coordinate 0 of max_pnt is greater than " << half_box_size;
     ASSERT_GE(half_box_size, max_pnt[1]) << "coordinate 1 of max_pnt is greater than " << half_box_size;
     ASSERT_GE(half_box_size, max_pnt[2]) << "coordinate 2 of max_pnt is greater than " << half_box_size;

     for (auto& point : pnts_list)
     {
         delete point;
     }
}

TEST(GeoLibAABB, RandomNumberOfPointsRandomPointInAList)
{
    /* initialize random seed: */
     srand ( static_cast<unsigned>(time(nullptr)) );
     int n (rand() % 1000000);
     int box_size (rand());
     double half_box_size(box_size/2);
     double minus_half_box_size(-half_box_size);

     // fill list with points
     std::list<GeoLib::Point> pnts_list;
     for (int k(0); k<n; k++) {
         pnts_list.emplace_back(rand() % box_size - half_box_size,
                                rand() % box_size - half_box_size,
                                rand() % box_size - half_box_size);
     }

     // construct from list points a axis algined bounding box
     GeoLib::AABB aabb(pnts_list.begin(), pnts_list.end());

     MathLib::Point3d const& min_pnt(aabb.getMinPoint());
     MathLib::Point3d const& max_pnt(aabb.getMaxPoint());

     ASSERT_LE(minus_half_box_size, min_pnt[0]) << "coordinate 0 of min_pnt is smaller than " << minus_half_box_size;
     ASSERT_LE(minus_half_box_size, min_pnt[1]) << "coordinate 1 of min_pnt is smaller than " << minus_half_box_size;
     ASSERT_LE(minus_half_box_size, min_pnt[2]) << "coordinate 2 of min_pnt is smaller than " << minus_half_box_size;

    // since the interval is half-open we have to move the upper bound also
     half_box_size = std::nexttoward(half_box_size, std::numeric_limits<double>::max());
     ASSERT_GE(half_box_size, max_pnt[0]) << "coordinate 0 of max_pnt is greater than " << half_box_size;
     ASSERT_GE(half_box_size, max_pnt[1]) << "coordinate 1 of max_pnt is greater than " << half_box_size;
     ASSERT_GE(half_box_size, max_pnt[2]) << "coordinate 2 of max_pnt is greater than " << half_box_size;
}

TEST(GeoLibAABB, RandomNumberOfPointersToRandomPointsInAVector)
{
    /* initialize random seed: */
     srand ( static_cast<unsigned>(time(nullptr)) );
     int n (rand() % 100000);
     int box_size (rand());
     double half_box_size(box_size/2);
     double minus_half_box_size(-half_box_size);

     // fill list with points
     std::vector<GeoLib::Point*> pnts;
     for (int k(0); k<n; k++) {
         pnts.push_back(new GeoLib::Point(rand() % box_size - half_box_size, rand() % box_size - half_box_size, rand() % box_size - half_box_size));
     }

     // construct from list points a axis aligned bounding box
     GeoLib::AABB aabb(pnts.begin(), pnts.end());

     MathLib::Point3d const& min_pnt(aabb.getMinPoint());
     MathLib::Point3d const& max_pnt(aabb.getMaxPoint());

     ASSERT_LE(minus_half_box_size, min_pnt[0]) << "coordinate 0 of min_pnt is smaller than " << minus_half_box_size;
     ASSERT_LE(minus_half_box_size, min_pnt[1]) << "coordinate 1 of min_pnt is smaller than " << minus_half_box_size;
     ASSERT_LE(minus_half_box_size, min_pnt[2]) << "coordinate 2 of min_pnt is smaller than " << minus_half_box_size;

    // since the interval is half-open we have to move the upper bound also
     half_box_size = std::nexttoward(half_box_size, std::numeric_limits<double>::max());
     ASSERT_GE(half_box_size, max_pnt[0]) << "coordinate 0 of max_pnt is greater than " << half_box_size;
     ASSERT_GE(half_box_size, max_pnt[1]) << "coordinate 1 of max_pnt is greater than " << half_box_size;
     ASSERT_GE(half_box_size, max_pnt[2]) << "coordinate 2 of max_pnt is greater than " << half_box_size;

     for (auto& point : pnts)
     {
         delete point;
     }
}

TEST(GeoLibAABB, RandomNumberOfPointsRandomPointInAVector)
{
    /* initialize random seed: */
     srand ( static_cast<unsigned>(time(nullptr)) );
     int n (rand() % 1000000);
     int box_size (rand());
     double half_box_size(box_size/2);
     double minus_half_box_size(-half_box_size);

     // fill list with points
     std::list<GeoLib::Point> pnts;
     for (int k(0); k<n; k++) {
         pnts.emplace_back(rand() % box_size - half_box_size,
                           rand() % box_size - half_box_size,
                           rand() % box_size - half_box_size);
     }

     // construct from list points a axis algined bounding box
     GeoLib::AABB aabb(pnts.begin(), pnts.end());

     MathLib::Point3d const& min_pnt(aabb.getMinPoint());
     MathLib::Point3d const& max_pnt(aabb.getMaxPoint());

     ASSERT_LE(minus_half_box_size, min_pnt[0]) << "coordinate 0 of min_pnt is smaller than " << minus_half_box_size;
     ASSERT_LE(minus_half_box_size, min_pnt[1]) << "coordinate 1 of min_pnt is smaller than " << minus_half_box_size;
     ASSERT_LE(minus_half_box_size, min_pnt[2]) << "coordinate 2 of min_pnt is smaller than " << minus_half_box_size;

    // since the interval is half-open we have to move the upper bound also
     half_box_size = std::nexttoward(half_box_size, std::numeric_limits<double>::max());
     ASSERT_GE(half_box_size, max_pnt[0]) << "coordinate 0 of max_pnt is greater than " << half_box_size;
     ASSERT_GE(half_box_size, max_pnt[1]) << "coordinate 1 of max_pnt is greater than " << half_box_size;
     ASSERT_GE(half_box_size, max_pnt[2]) << "coordinate 2 of max_pnt is greater than " << half_box_size;
}

TEST(GeoLibAABB, RandomNumberOfPointsRandomBox)
{
    /* initialize random seed: */
     srand (static_cast<unsigned>(time(nullptr)));
     int n (rand() % 1000000);
     int box_size_x (rand());
     int box_size_y (rand());
     int box_size_z (rand());
     double half_box_size_x(box_size_x/2);
     double half_box_size_y(box_size_y/2);
     double half_box_size_z(box_size_z/2);
     auto minus_half_box_size_x(static_cast<int>(-half_box_size_x));
     auto minus_half_box_size_y(static_cast<int>(-half_box_size_y));
     auto minus_half_box_size_z(static_cast<int>(-half_box_size_z));

     // fill list with points
     std::list<GeoLib::Point> pnts;
     for (int k(0); k<n; k++) {
         pnts.emplace_back(rand() % box_size_x - half_box_size_x,
                           rand() % box_size_y - half_box_size_y,
                           rand() % box_size_z - half_box_size_z);
     }

     // construct from list points a axis aligned bounding box
     GeoLib::AABB aabb(pnts.begin(), pnts.end());

     MathLib::Point3d const& min_pnt(aabb.getMinPoint());
     MathLib::Point3d const& max_pnt(aabb.getMaxPoint());

     ASSERT_LE(minus_half_box_size_x, min_pnt[0]) << "coordinate 0 of min_pnt is smaller than " << minus_half_box_size_x;
     ASSERT_LE(minus_half_box_size_y, min_pnt[1]) << "coordinate 1 of min_pnt is smaller than " << minus_half_box_size_y;
     ASSERT_LE(minus_half_box_size_z, min_pnt[2]) << "coordinate 2 of min_pnt is smaller than " << minus_half_box_size_z;

    // since the interval is half-open we have to move the upper bound also
     half_box_size_x = std::nexttoward(half_box_size_x, std::numeric_limits<double>::max());
     half_box_size_y = std::nexttoward(half_box_size_y, std::numeric_limits<double>::max());
     half_box_size_z = std::nexttoward(half_box_size_z, std::numeric_limits<double>::max());
     ASSERT_GE(half_box_size_x, max_pnt[0]) << "coordinate 0 of max_pnt is greater than " << half_box_size_x;
     ASSERT_GE(half_box_size_y, max_pnt[1]) << "coordinate 1 of max_pnt is greater than " << half_box_size_y;
     ASSERT_GE(half_box_size_z, max_pnt[2]) << "coordinate 2 of max_pnt is greater than " << half_box_size_z;
}

TEST(GeoLib, AABBAllPointsWithNegativeCoordinatesI)
{
    std::vector<GeoLib::Point*> pnts;
    pnts.push_back(new GeoLib::Point(-1, -1, -1));
    pnts.push_back(new GeoLib::Point(-10, -10, -10));

    std::vector<std::size_t> ids;
    ids.push_back(0);
    ids.push_back(1);
    GeoLib::AABB aabb(pnts, ids);

    MathLib::Point3d const& max_pnt(aabb.getMaxPoint());

    ASSERT_NEAR(-1.0, max_pnt[0], std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(-1.0, max_pnt[1], std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(-1.0, max_pnt[2], std::numeric_limits<double>::epsilon());

    for (auto p : pnts)
        delete p;
}

TEST(GeoLib, AABBAllPointsWithNegativeCoordinatesII)
{
    std::vector<GeoLib::Point> pnts;

    pnts.emplace_back(-1, -1, -1);
    pnts.emplace_back(-10, -10, -10);

    // construct from points of the vector a axis aligned bounding box
    GeoLib::AABB aabb(pnts.begin(), pnts.end());

    MathLib::Point3d const& max_pnt(aabb.getMaxPoint());

    ASSERT_NEAR(-1.0, max_pnt[0], std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(-1.0, max_pnt[1], std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(-1.0, max_pnt[2], std::numeric_limits<double>::epsilon());
}

// This test creates an AABB containing a single point.
// This point is moved in every space direction to the next possible double value.
// It is checked that the AABB contains the point iff it has not been moved.
TEST(GeoLib, AABBSinglePoint)
{
    std::random_device rd;
    std::mt19937 random_number_generator(rd());
    std::uniform_real_distribution<double> rnd(-1000000.0, 10000000.0);
    std::vector<GeoLib::Point> pnts;
    GeoLib::Point p{{{rnd(random_number_generator),
        rnd(random_number_generator),rnd(random_number_generator)}}};
    pnts.push_back(p);

    ASSERT_EQ(1u, pnts.size());

    // construct from points of the vector a axis aligned bounding box
    GeoLib::AABB aabb(pnts.begin(), pnts.end());

    double const to_lowest(std::numeric_limits<double>::lowest());
    double const to_max(std::numeric_limits<double>::max());

    // Check the point within the aabb (i==j==k). The outer 26 (3 * 3 * 3 - 1)
    // points around the aabb are also checked.
    for (int i(-1); i<2; ++i) {
        // Modify the first coordinate of p.
        if (i==-1) p[0] = std::nextafter(pnts.front()[0], to_lowest);
        else if (i==0) p[0] = pnts.front()[0];
        else p[0] = std::nextafter(pnts.front()[0], to_max);
        for (int j(-1); j<2; ++j) {
            // Modify the second coordinate of p.
            if (j==-1) p[1] = std::nextafter(pnts.front()[1], to_lowest);
            else if (j==0) p[1] = pnts.front()[1];
            else p[1] = std::nextafter(pnts.front()[1], to_max);
            for (int k(-1); k<2; ++k) {
                // Modify the third coordinate of p.
                if (k==-1) p[2] = std::nextafter(pnts.front()[2], to_lowest);
                else if (k==0) p[2] = pnts.front()[2];
                else p[2] = std::nextafter(pnts.front()[2], to_max);
                if (i == 0 && j == 0 && k == 0)
                    ASSERT_TRUE(aabb.containsPoint(p, 0));
                else
                    ASSERT_FALSE(aabb.containsPoint(p, 0));
            }
        }
    }
}

