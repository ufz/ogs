/**
 * \file
 * \brief  Tests of Raster class.
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gtest/gtest.h>

#include "GeoLib/Raster.h"

TEST(GeoLibRaster, EmptyRaster)
{
    GeoLib::RasterHeader const header{0, 0,    0, MathLib::Point3d{{0, 0, 0}},
                                      0, -9999};
    std::vector<double> const data;
    GeoLib::Raster const empty_raster{header, data.begin(), data.end()};
    // accessing any pixel in the raster should throw an exception
    EXPECT_THROW(empty_raster(0, 0), std::runtime_error);
}

TEST(GeoLibRaster, OnePixelRasterEmptyData)
{
    GeoLib::RasterHeader const header{1, 1,    0, MathLib::Point3d{{0, 0, 0}},
                                      0, -9999};
    std::vector<double> const empty_data;
    EXPECT_THROW(GeoLib::Raster(header, empty_data.begin(), empty_data.end()),
                 std::out_of_range);
}

TEST(GeoLibRaster, OnePixelRaster)
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

TEST(GeoLibRaster, OneRowByTwoColumnsRaster)
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

TEST(GeoLibRaster, TwoRowsByFourColumnsRaster)
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
