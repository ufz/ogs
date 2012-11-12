/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 * \file TestAABB.cpp
 *
 *  Created on  Oct 24, 2012 by Thomas Fischer
 */

// ** INCLUDES **
#include "gtest/gtest.h"

// std lib
#include <cstdlib>
#include <ctime>
#include <list>
#include <iostream>

#include "Point.h"
#include "AABB.h"

TEST(GeoLibAABB, RandomNumberOfPointsPointersToRandomPointsInAList)
{
	/* initialize random seed: */
	 srand ( time(NULL) );
	 int n (rand() % 100000);
	 int box_size (rand());
	 int half_box_size(box_size/2);
	 int minus_half_box_size(-1.0 * half_box_size);

	 // fill list with points
	 std::list<GeoLib::Point*> pnts_list;
	 for (int k(0); k<n; k++) {
		 pnts_list.push_back(new GeoLib::Point(rand() % box_size - half_box_size, rand() % box_size - half_box_size, rand() % box_size - half_box_size));
	 }

	 std::cout << "testing with " << n << " points" << std::endl;

	 // construct from list points a axis algined bounding box
	 GeoLib::AABB<GeoLib::Point> aabb(pnts_list.begin(), pnts_list.end());

	 GeoLib::Point const& min_pnt(aabb.getMinPoint());
	 GeoLib::Point const& max_pnt(aabb.getMaxPoint());

	 ASSERT_LE(minus_half_box_size, min_pnt[0]) << "coordinate 0 of min_pnt is smaller than " << minus_half_box_size;
	 ASSERT_LE(minus_half_box_size, min_pnt[1]) << "coordinate 1 of min_pnt is smaller than " << minus_half_box_size;
	 ASSERT_LE(minus_half_box_size, min_pnt[2]) << "coordinate 2 of min_pnt is smaller than " << minus_half_box_size;
	 ASSERT_GE(half_box_size, max_pnt[0]) << "coordinate 0 of max_pnt is greater than " << half_box_size;
	 ASSERT_GE(half_box_size, max_pnt[1]) << "coordinate 1 of max_pnt is greater than " << half_box_size;
	 ASSERT_GE(half_box_size, max_pnt[2]) << "coordinate 2 of max_pnt is greater than " << half_box_size;

	 for (std::list<GeoLib::Point*>::iterator it(pnts_list.begin()); it != pnts_list.end(); it++) {
		 delete *it;
	 }
}

TEST(GeoLibAABB, RandomNumberOfPointsRandomPointInAList)
{
	/* initialize random seed: */
	 srand ( time(NULL) );
	 int n (rand() % 1000000);
	 int box_size (rand());
	 int half_box_size(box_size/2);
	 int minus_half_box_size(-1.0 * half_box_size);

	 // fill list with points
	 std::list<GeoLib::Point> pnts_list;
	 for (int k(0); k<n; k++) {
		 pnts_list.push_back(GeoLib::Point(rand() % box_size - half_box_size, rand() % box_size - half_box_size, rand() % box_size - half_box_size));
	 }

	 std::cout << "testing with " << n << " points" << std::endl;

	 // construct from list points a axis algined bounding box
	 GeoLib::AABB<GeoLib::Point> aabb(pnts_list.begin(), pnts_list.end());

	 GeoLib::Point const& min_pnt(aabb.getMinPoint());
	 GeoLib::Point const& max_pnt(aabb.getMaxPoint());

	 ASSERT_LE(minus_half_box_size, min_pnt[0]) << "coordinate 0 of min_pnt is smaller than " << minus_half_box_size;
	 ASSERT_LE(minus_half_box_size, min_pnt[1]) << "coordinate 1 of min_pnt is smaller than " << minus_half_box_size;
	 ASSERT_LE(minus_half_box_size, min_pnt[2]) << "coordinate 2 of min_pnt is smaller than " << minus_half_box_size;
	 ASSERT_GE(half_box_size, max_pnt[0]) << "coordinate 0 of max_pnt is greater than " << half_box_size;
	 ASSERT_GE(half_box_size, max_pnt[1]) << "coordinate 1 of max_pnt is greater than " << half_box_size;
	 ASSERT_GE(half_box_size, max_pnt[2]) << "coordinate 2 of max_pnt is greater than " << half_box_size;
}

TEST(GeoLibAABB, RandomNumberOfPointsPointersToRandomPointsInAVector)
{
	/* initialize random seed: */
	 srand ( time(NULL) );
	 int n (rand() % 100000);
	 int box_size (rand());
	 int half_box_size(box_size/2);
	 int minus_half_box_size(-1.0 * half_box_size);

	 // fill list with points
	 std::vector<GeoLib::Point*> pnts;
	 for (int k(0); k<n; k++) {
		 pnts.push_back(new GeoLib::Point(rand() % box_size - half_box_size, rand() % box_size - half_box_size, rand() % box_size - half_box_size));
	 }

	 std::cout << "testing with " << n << " points" << std::endl;

	 // construct from list points a axis algined bounding box
	 GeoLib::AABB<GeoLib::Point> aabb(pnts.begin(), pnts.end());

	 GeoLib::Point const& min_pnt(aabb.getMinPoint());
	 GeoLib::Point const& max_pnt(aabb.getMaxPoint());

	 ASSERT_LE(minus_half_box_size, min_pnt[0]) << "coordinate 0 of min_pnt is smaller than " << minus_half_box_size;
	 ASSERT_LE(minus_half_box_size, min_pnt[1]) << "coordinate 1 of min_pnt is smaller than " << minus_half_box_size;
	 ASSERT_LE(minus_half_box_size, min_pnt[2]) << "coordinate 2 of min_pnt is smaller than " << minus_half_box_size;
	 ASSERT_GE(half_box_size, max_pnt[0]) << "coordinate 0 of max_pnt is greater than " << half_box_size;
	 ASSERT_GE(half_box_size, max_pnt[1]) << "coordinate 1 of max_pnt is greater than " << half_box_size;
	 ASSERT_GE(half_box_size, max_pnt[2]) << "coordinate 2 of max_pnt is greater than " << half_box_size;

	 for (std::vector<GeoLib::Point*>::iterator it(pnts.begin()); it != pnts.end(); it++) {
		 delete *it;
	 }
}

TEST(GeoLibAABB, RandomNumberOfPointsRandomPointInAVector)
{
	/* initialize random seed: */
	 srand ( time(NULL) );
	 int n (rand() % 1000000);
	 int box_size (rand());
	 int half_box_size(box_size/2);
	 int minus_half_box_size(-1.0 * half_box_size);

	 // fill list with points
	 std::list<GeoLib::Point> pnts;
	 for (int k(0); k<n; k++) {
		 pnts.push_back(GeoLib::Point(rand() % box_size - half_box_size, rand() % box_size - half_box_size, rand() % box_size - half_box_size));
	 }

	 std::cout << "testing with " << n << " points" << std::endl;

	 // construct from list points a axis algined bounding box
	 GeoLib::AABB<GeoLib::Point> aabb(pnts.begin(), pnts.end());

	 GeoLib::Point const& min_pnt(aabb.getMinPoint());
	 GeoLib::Point const& max_pnt(aabb.getMaxPoint());

	 ASSERT_LE(minus_half_box_size, min_pnt[0]) << "coordinate 0 of min_pnt is smaller than " << minus_half_box_size;
	 ASSERT_LE(minus_half_box_size, min_pnt[1]) << "coordinate 1 of min_pnt is smaller than " << minus_half_box_size;
	 ASSERT_LE(minus_half_box_size, min_pnt[2]) << "coordinate 2 of min_pnt is smaller than " << minus_half_box_size;
	 ASSERT_GE(half_box_size, max_pnt[0]) << "coordinate 0 of max_pnt is greater than " << half_box_size;
	 ASSERT_GE(half_box_size, max_pnt[1]) << "coordinate 1 of max_pnt is greater than " << half_box_size;
	 ASSERT_GE(half_box_size, max_pnt[2]) << "coordinate 2 of max_pnt is greater than " << half_box_size;
}

