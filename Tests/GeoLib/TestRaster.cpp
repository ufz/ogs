/**
 * \file
 * \brief  Tests of Raster class.
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gtest/gtest.h>

#include <numeric>
#include <random>

#include "GeoLib/Point.h"
#include "GeoLib/Raster.h"

TEST(GeoLibFixedRaster, EmptyRaster)
{
    GeoLib::RasterHeader const header{0, 0,    0, MathLib::Point3d{{0, 0, 0}},
                                      0, -9999};
    std::vector<double> const data;
    GeoLib::Raster const empty_raster{header, data.begin(), data.end()};
    // accessing any pixel in the raster should throw an exception
    EXPECT_THROW(empty_raster(0, 0), std::runtime_error);
}

TEST(GeoLibFixedRaster, OnePixelRasterEmptyData)
{
    GeoLib::RasterHeader const header{1, 1,    0, MathLib::Point3d{{0, 0, 0}},
                                      0, -9999};
    std::vector<double> const empty_data;
    EXPECT_THROW(GeoLib::Raster(header, empty_data.begin(), empty_data.end()),
                 std::out_of_range);
}

TEST(GeoLibFixedRaster, OnePixelRaster)
{
    GeoLib::RasterHeader const header{1, 1,    0, MathLib::Point3d{{0, 0, 0}},
                                      0, -9999};
    std::vector<double> const data{1};
    GeoLib::Raster const one_pixel_raster{header, data.begin(), data.end()};
    EXPECT_EQ(1.0, one_pixel_raster(0, 0));
    // accessing any other pixel in the raster should throw an exception
    EXPECT_THROW(one_pixel_raster(1, 0), std::runtime_error);
    EXPECT_THROW(one_pixel_raster(0, 1), std::runtime_error);
    EXPECT_THROW(one_pixel_raster(1, 1), std::runtime_error);
    EXPECT_THROW(one_pixel_raster(-1, 0), std::runtime_error);
}

TEST(GeoLibFixedRaster, OneRowByTwoColumnsRaster)
{
    GeoLib::RasterHeader const header{2, 1,    0, MathLib::Point3d{{0, 0, 0}},
                                      0, -9999};
    std::vector<double> const data{1, 2};
    GeoLib::Raster const one_pixel_raster{header, data.begin(), data.end()};
    EXPECT_EQ(1.0, one_pixel_raster(0, 0));
    EXPECT_EQ(2.0, one_pixel_raster(0, 1));
    // accessing any other pixel in the raster should throw an exception
    EXPECT_THROW(one_pixel_raster(1, 0), std::runtime_error);
    EXPECT_THROW(one_pixel_raster(1, 1), std::runtime_error);
}

TEST(GeoLibFixedRaster, TwoRowsByFourColumnsRaster)
{
    GeoLib::RasterHeader const header{4, 2,    0, MathLib::Point3d{{0, 0, 0}},
                                      0, -9999};
    std::vector<double> const data{1, 2, 3, 4, 5, 6, 7, 8};
    GeoLib::Raster const raster{header, data.begin(), data.end()};
    // compare first row
    for (std::size_t c = 0; c < header.n_cols; ++c)
    {
        EXPECT_EQ(double(c + 1 + header.n_cols), raster(0, c));
    }
    // compare second row
    for (std::size_t c = 0; c < header.n_cols; ++c)
    {
        EXPECT_EQ(double(c + 1), raster(1, c));
    }
}

struct GeoLibRaster : public testing::Test
{
    GeoLibRaster()
    {
        std::random_device random_device;
        std::default_random_engine random_engine(random_device());
        std::uniform_int_distribution<std::size_t> uniform_dist(1, 20);

        n_columns = uniform_dist(random_engine);
        n_rows = uniform_dist(random_engine);

        std::normal_distribution<> normal_dist(100, 90);
        cell_size = normal_dist(random_engine);
        header = GeoLib::RasterHeader{n_columns, n_rows,
                                      0,         MathLib::Point3d{{0, 0, 0}},
                                      cell_size, -9999};
        data.resize(n_columns * n_rows);
        std::iota(data.begin(), data.end(), 1);
    }

    GeoLib::RasterHeader header;
    std::size_t n_columns;
    std::size_t n_rows;
    double cell_size;
    std::vector<double> data;
};

TEST_F(GeoLibRaster, CopyConstructor)
{
    GeoLib::Raster const raster{header, data.begin(), data.end()};
    GeoLib::Raster const raster_copy(raster);

    ASSERT_EQ(raster, raster_copy);
}

TEST_F(GeoLibRaster, AssignmentOperator)
{
    GeoLib::Raster const raster{header, data.begin(), data.end()};

    std::vector<double> const empty_data{};
    GeoLib::RasterHeader const dummy_header{
        0, 0, 0, GeoLib::Point{{0.0, 0.0, 0.0}}, 0, 0};

    GeoLib::Raster raster_copy{dummy_header, empty_data.begin(),
                               empty_data.end()};
    raster_copy = raster;
    ASSERT_EQ(raster, raster_copy);
}

TEST_F(GeoLibRaster, MoveAssignmentOperator)
{
    GeoLib::Raster raster{header, data.begin(), data.end()};
    GeoLib::Raster const raster_copy{header, data.begin(), data.end()};

    GeoLib::Raster raster_move = std::move(raster);
    ASSERT_EQ(raster_copy, raster_move);
}

TEST_F(GeoLibRaster, MoveConstructor)
{
    GeoLib::Raster raster{header, data.begin(), data.end()};
    GeoLib::Raster const raster_copy{header, data.begin(), data.end()};

    GeoLib::Raster raster_move{std::move(raster)};
    ASSERT_EQ(raster_copy, raster_move);
}
