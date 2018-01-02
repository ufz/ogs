/**
 * @file TestGrid.cpp
 * @author Thomas Fischer
 * @date Mar 22, 2013
 * @brief Tests for class GeoLib::Grid.
 *
 * @copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <memory>

#include "gtest/gtest.h"

#include "GeoLib/Point.h"
#include "GeoLib/Grid.h"

#include "MathLib/MathTools.h"

TEST(GeoLib, InsertZeroPointsInGrid)
{
    std::vector<GeoLib::Point*> pnts;
    ASSERT_THROW(GeoLib::Grid<GeoLib::Point> grid(pnts.begin(), pnts.end()), std::runtime_error);
}

TEST(GeoLib, InsertOnePointInGrid)
{
    std::vector<GeoLib::Point*> pnts;
    pnts.push_back(new GeoLib::Point(0.0,0.0,0.0));
    ASSERT_NO_THROW(GeoLib::Grid<GeoLib::Point> grid(pnts.begin(), pnts.end()));
    std::for_each(pnts.begin(), pnts.end(), std::default_delete<GeoLib::Point>());
}

TEST(GeoLib, InsertTwoPointsInGrid)
{
    std::vector<GeoLib::Point*> pnts;
    pnts.push_back(new GeoLib::Point(4.5, -400.0, 0.0));
    pnts.push_back(new GeoLib::Point(50, -300.0, 0.0));
    ASSERT_NO_THROW(GeoLib::Grid<GeoLib::Point> grid(pnts.begin(), pnts.end()));
    std::for_each(pnts.begin(), pnts.end(), std::default_delete<GeoLib::Point>());
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
    std::for_each(pnts.begin(), pnts.end(), std::default_delete<GeoLib::Point>());
}

TEST(GeoLib, InsertPointsWithSameXCoordinateInGrid)
{
    std::vector<GeoLib::Point*> pnts;
    pnts.push_back(new GeoLib::Point(0,0,0));
    pnts.push_back(new GeoLib::Point(0,1,0));
    pnts.push_back(new GeoLib::Point(0,1,0.1));

    GeoLib::Grid<GeoLib::Point> grid(pnts.begin(), pnts.end());
    std::for_each(pnts.begin(), pnts.end(), std::default_delete<GeoLib::Point>());
}

TEST(GeoLib, SearchNearestPointInGrid)
{
    std::vector<GeoLib::Point*> pnts;
    pnts.push_back(new GeoLib::Point(0.0,0.0,0.0));
    GeoLib::Grid<GeoLib::Point> *grid(nullptr);
    ASSERT_NO_THROW(grid = new GeoLib::Grid<GeoLib::Point>(pnts.begin(), pnts.end()));

    GeoLib::Point p0(0,10,10);
    GeoLib::Point* res(grid->getNearestPoint(p0));
    ASSERT_EQ(0.0, sqrt(MathLib::sqrDist(*res, *pnts[0])));

    delete grid;
    std::for_each(pnts.begin(), pnts.end(), std::default_delete<GeoLib::Point>());
}

TEST(GeoLib, SearchNearestPointsInDenseGrid)
{
    const std::size_t i_max(50), j_max(50), k_max(50);
    std::vector<GeoLib::Point*> pnts(i_max*j_max*k_max);

    // fill the vector with equi-distant points in the
    // cube [0,(i_max-1)/i_max] x [0,(j_max-1)/j_max] x [0,(k_max-1)/k_max]
    for (std::size_t i(0); i < i_max; i++) {
        std::size_t offset0(i * j_max * k_max);
        for (std::size_t j(0); j < j_max; j++) {
            std::size_t offset1(j * k_max + offset0);
            for (std::size_t k(0); k < k_max; k++) {
                pnts[offset1 + k] = new GeoLib::Point(
                    std::array<double,3>({{static_cast<double>(i) / i_max,
                        static_cast<double>(j) / j_max,
                        static_cast<double>(k) / k_max}}), offset1+k);
            }
        }
    }

    // create the grid
    GeoLib::Grid<GeoLib::Point>* grid(nullptr);
    ASSERT_NO_THROW(grid = new GeoLib::Grid<GeoLib::Point> (pnts.begin(), pnts.end()));

    // search point (1,1,1) is outside of the point set
    GeoLib::Point search_pnt(std::array<double,3>({{1,1,1}}), 0);
    GeoLib::Point* res(grid->getNearestPoint(search_pnt));
    ASSERT_EQ(i_max*j_max*k_max-1, res->getID());
    ASSERT_NEAR(sqrt(3.0)/50.0, sqrt(MathLib::sqrDist(*res, search_pnt)), std::numeric_limits<double>::epsilon());

    // search point (0,1,1) is outside of the point set
    search_pnt[0] = 0;
    res = grid->getNearestPoint(search_pnt);
    ASSERT_EQ(j_max*k_max - 1, res->getID());
    ASSERT_NEAR(sqrt(2.0)/50.0, sqrt(MathLib::sqrDist(*res, search_pnt)), std::numeric_limits<double>::epsilon());

    // search point (0.5,1,1) is outside of the point set
    search_pnt[0] = 0.5;
    res = grid->getNearestPoint(search_pnt);
    ASSERT_EQ(j_max*k_max*(i_max/2 + 1) - 1, res->getID());
    ASSERT_NEAR(sqrt(2.0)/50.0, sqrt(MathLib::sqrDist(*res, search_pnt)), std::numeric_limits<double>::epsilon());

    // checking only every fourth point per direction to reduce the run time of
    // the test
    for (std::size_t i(0); i < i_max; i=i+4) {
        std::size_t offset0(i * j_max * k_max);
        for (std::size_t j(0); j < j_max; j=j+4) {
            std::size_t offset1(j * k_max + offset0);
            for (std::size_t k(0); k < k_max; k=k+4) {
                res = grid->getNearestPoint(*pnts[offset1+k]);
                ASSERT_EQ(offset1+k, res->getID());
                ASSERT_NEAR(sqrt(MathLib::sqrDist(*res, *pnts[offset1+k])), 0.0, std::numeric_limits<double>::epsilon());
            }
        }
    }

    delete grid;
    std::for_each(pnts.begin(), pnts.end(), std::default_delete<GeoLib::Point>());
}
