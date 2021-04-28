/**
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 */

#include <gtest/gtest.h>
#include <vtkImageAlgorithm.h>
#include <vtkImageData.h>
#include <vtkPNGReader.h>
#include <vtkSmartPointer.h>

#include <memory>

#include "Applications/DataExplorer/VtkVis/VtkRaster.h"
#include "Applications/FileIO/AsciiRasterInterface.h"
#include "GeoLib/Raster.h"
#include "InfoLib/TestInfo.h"

TEST(TestVtkRaster, TestPNGReader)
{
    std::string name =
        TestInfoLib::TestInfo::data_path + "/FileIO/testraster.png";
    vtkSmartPointer<vtkImageAlgorithm> img = VtkRaster::loadImage(name);
    img->Update();
    EXPECT_TRUE(img != nullptr);
    EXPECT_TRUE(img->GetOutput() != nullptr);
    double val[3];
    img->GetOutput()->GetSpacing(val);
    for (double i : val)
    {
        EXPECT_NEAR(60, i, std::numeric_limits<double>::epsilon());
    }
    img->GetOutput()->GetOrigin(val);
    EXPECT_NEAR(5000, val[0], std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(-400, val[1], std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(0, val[2], std::numeric_limits<double>::epsilon());
    int extent[3];
    img->GetOutput()->GetDimensions(extent);
    EXPECT_EQ(12, extent[0]);
    EXPECT_EQ(10, extent[1]);
    EXPECT_EQ(1, extent[2]);
}

TEST(TestVtkRaster, TestASCReader)
{
    std::string name = TestInfoLib::TestInfo::data_path +
                       "/MeshGeoToolsLib/Hamburg/00-raster.asc";
    vtkSmartPointer<vtkImageAlgorithm> img = VtkRaster::loadImage(name);
    img->Update();
    EXPECT_TRUE(img != nullptr);
    EXPECT_TRUE(img->GetOutput() != nullptr);
    std::unique_ptr<GeoLib::Raster> raster(
        FileIO::AsciiRasterInterface::getRasterFromASCFile(name));
    EXPECT_EQ(298, raster->getHeader().n_cols);
    EXPECT_EQ(205, raster->getHeader().n_rows);
    EXPECT_NEAR(25, raster->getHeader().cell_size,
                std::numeric_limits<double>::epsilon());

    double val[3];
    img->GetOutput()->GetSpacing(val);
    for (double i : val)
    {
        EXPECT_NEAR(raster->getHeader().cell_size, i,
                    std::numeric_limits<double>::epsilon());
    }

    double const half_pix = raster->getHeader().cell_size / 2.0;
    img->GetOutput()->GetOrigin(val);
    EXPECT_NEAR(raster->getHeader().origin[0], val[0] - half_pix,
                std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(raster->getHeader().origin[1], val[1] - half_pix,
                std::numeric_limits<double>::epsilon());
    EXPECT_NEAR(raster->getHeader().origin[2], val[2],
                std::numeric_limits<double>::epsilon());

    int extent[3];
    img->GetOutput()->GetDimensions(extent);
    EXPECT_EQ(raster->getHeader().n_cols, extent[0]);
    EXPECT_EQ(raster->getHeader().n_rows, extent[1]);
    EXPECT_EQ(raster->getHeader().n_depth, extent[2]);
}
