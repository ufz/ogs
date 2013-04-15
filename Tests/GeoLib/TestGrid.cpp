/**
 * @file TestGrid.cpp
 * @author Thomas Fischer
 * @date Mar 22, 2013
 * @brief 
 *
 * @copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include "gtest/gtest.h"

#include "GeoLib/Point.h"
#include "GeoLib/PointWithID.h"
#include "GeoLib/Grid.h"

#include "MathTools.h"

TEST(GeoLib, InsertZeroPointsInGrid)
{
	std::vector<GeoLib::Point*> pnts;
	ASSERT_THROW(GeoLib::Grid<GeoLib::Point> grid(pnts.begin(), pnts.end()), std::invalid_argument);
}

TEST(GeoLib, InsertOnePointInGrid)
{
	std::vector<GeoLib::Point*> pnts;
	pnts.push_back(new GeoLib::Point(0.0,0.0,0.0));
	ASSERT_NO_THROW(GeoLib::Grid<GeoLib::Point> grid(pnts.begin(), pnts.end()));
}

TEST(GeoLib, InsertManyPointsInGrid)
{
	const std::size_t i_max(100), j_max(100), k_max(100);
	std::vector<GeoLib::Point*> pnts(i_max*j_max*k_max);

	// fill the vector with points
	for (std::size_t i(0); i < i_max; i++) {
		std::size_t offset0(i * j_max * k_max);
		for (std::size_t j(0); j < j_max; j++) {
			std::size_t offset1(j * k_max + offset0);
			for (std::size_t k(0); k < k_max; k++) {
				pnts[offset1 + k] = new GeoLib::Point(static_cast<double>(i) / i_max,
						static_cast<double>(j) / j_max, static_cast<double>(k) / k_max);
			}
		}
	}

	ASSERT_NO_THROW(GeoLib::Grid<GeoLib::Point> grid(pnts.begin(), pnts.end()));
}

TEST(GeoLib, SearchNearestPointInGrid)
{
	std::vector<GeoLib::Point*> pnts;
	pnts.push_back(new GeoLib::Point(0.0,0.0,0.0));
	GeoLib::Grid<GeoLib::Point> *grid(nullptr);
	ASSERT_NO_THROW(grid = new GeoLib::Grid<GeoLib::Point>(pnts.begin(), pnts.end()));

	GeoLib::Point p0(0,10,10);
	GeoLib::Point* res(grid->getNearestPoint(p0.getCoords()));
	ASSERT_EQ(sqrt(MathLib::sqrDist(res, pnts[0])), 0.0);
}

TEST(GeoLib, SearchNearestPointsInDenseGrid)
{
	const std::size_t i_max(50), j_max(50), k_max(50);
	std::vector<GeoLib::PointWithID*> pnts(i_max*j_max*k_max);

	// fill the vector with equi-distant points in the
	// cube [0,(i_max-1)/i_max] x [0,(j_max-1)/j_max] x [0,(k_max-1)/k_max]
	for (std::size_t i(0); i < i_max; i++) {
		std::size_t offset0(i * j_max * k_max);
		for (std::size_t j(0); j < j_max; j++) {
			std::size_t offset1(j * k_max + offset0);
			for (std::size_t k(0); k < k_max; k++) {
				pnts[offset1 + k] = new GeoLib::PointWithID(static_cast<double>(i) / i_max,
						static_cast<double>(j) / j_max, static_cast<double>(k) / k_max, offset1+k);
			}
		}
	}

	// create the grid
	GeoLib::Grid<GeoLib::PointWithID>* grid(nullptr);
	ASSERT_NO_THROW(grid = new GeoLib::Grid<GeoLib::PointWithID> (pnts.begin(), pnts.end()));

	// search point (1,1,1) is outside of the point set
	GeoLib::PointWithID search_pnt(1,1,1, 0);
	GeoLib::PointWithID* res(grid->getNearestPoint(search_pnt.getCoords()));
	ASSERT_EQ(res->getID(), i_max*j_max*k_max-1);
	ASSERT_NEAR(sqrt(MathLib::sqrDist(res, &search_pnt)), sqrt(3.0)/50.0, std::numeric_limits<double>::epsilon());

	// search point (0,1,1) is outside of the point set
	search_pnt[0] = 0;
	res = grid->getNearestPoint(search_pnt.getCoords());
	ASSERT_EQ(res->getID(), j_max*k_max - 1);
	ASSERT_NEAR(sqrt(MathLib::sqrDist(res, &search_pnt)), sqrt(2.0)/50.0, std::numeric_limits<double>::epsilon());

	// search point (0.5,1,1) is outside of the point set
	search_pnt[0] = 0.5;
	res = grid->getNearestPoint(search_pnt.getCoords());
	ASSERT_EQ(res->getID(), j_max*k_max*(i_max/2 + 1) - 1);
	ASSERT_NEAR(sqrt(MathLib::sqrDist(res, &search_pnt)), sqrt(2.0)/50.0, std::numeric_limits<double>::epsilon());

	for (std::size_t i(0); i < i_max; i++) {
		std::size_t offset0(i * j_max * k_max);
		for (std::size_t j(0); j < j_max; j++) {
			std::size_t offset1(j * k_max + offset0);
			for (std::size_t k(0); k < k_max; k++) {
				res = grid->getNearestPoint(pnts[offset1+k]->getCoords());
				ASSERT_EQ(res->getID(), offset1+k);
				ASSERT_NEAR(sqrt(MathLib::sqrDist(res, pnts[offset1+k])), 0.0, std::numeric_limits<double>::min());
			}
		}
	}

}
