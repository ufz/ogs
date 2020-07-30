/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <algorithm>
#include <memory>
#include <string>

// ThirdParty
#include <tclap/CmdLine.h>

#include "InfoLib/GitInfo.h"

#include <vtkPolyDataAlgorithm.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkSmartPointer.h>

#include "BaseLib/FileTools.h"

#include "Applications/FileIO/AsciiRasterInterface.h"
#include "Applications/DataExplorer/VtkVis/VtkGeoImageSource.h"
#include "Applications/DataExplorer/VtkVis/VtkImageDataToPointCloudFilter.h"

int main(int argc, char *argv[])
{
    TCLAP::CmdLine cmd(
        "Batch conversion of raster files into VTK point cloud data.\n"
        "A series of raster files is converted into point clouds by creating a "
        "number of random points based on the intensity value in each pixel "
        "(the higher the value, the more points. The [x,y]-extend of these "
        "random points is limited by the extent of the pixel they represent, "
        "the elevation is limited by the max-elevation parameter. Likewise, the "
        "maximum number of points per pixel is limited by the max-points "
        "parameter. The increase of the number of points in relation to the "
        "pixel value is indicated by gamma, where 1 is a linear increase and "
        "values smaller or larger than 1 indicate a logarithmic or exponential "
        "increase, respectively.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2020, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);
    TCLAP::ValueArg<double> gamma_arg(
        "g", "gamma",
        "Exponental increase coefficient. A value of '1' indicates a linear "
        "increase, smaller or larger values indicate a logarithmic or "
        "exponential increase, respectively.",
        false, 0.0, "exponential increase (default: 1)");
    cmd.add(gamma_arg);
    TCLAP::ValueArg<std::size_t> max_points_arg(
        "p", "max-points", "Maximum number of points per pixel.", false, 0,
        "maximum points per pixel (default: 1000)");
    cmd.add(max_points_arg);
    TCLAP::ValueArg<double> max_height_arg(
        "e", "max-elevation", "Maximum elevation of randomly created points.",
        false, 0.0, "maximum point elevation (default: 5000");
    cmd.add(max_height_arg);
    TCLAP::ValueArg<std::string> output_arg(
        "o", "output",
        "Path to output polydata files (*.vtp). Specifying 'file.vtp' will "
        "make the programm write a series of files called 'file0.vtp', "
        "'file1.vtp', etc.",
        true, "", "output file name");
    cmd.add(output_arg);
    TCLAP::ValueArg<std::string> input_arg(
        "i", "input",
        "Path to input raster files (*.asc). Specifying 'file.asc' will make "
        "the programm look for a series of files called 'file0.asc', "
        "'file1.asc', etc.",
        true, "", "input file name");
    cmd.add(input_arg);
    cmd.parse( argc, argv );

    std::string const input_name = input_arg.getValue().c_str();
    std::string const output_name = output_arg.getValue().c_str();

    double const max_height =
        (max_height_arg.isSet()) ? max_height_arg.getValue() : 5000;
    std::size_t const max_points =
        (max_points_arg.isSet()) ? max_points_arg.getValue() : 1000;
    double const gamma =
        (gamma_arg.isSet()) ? gamma_arg.getValue() : 1;

    std::string const ibase_name = BaseLib::dropFileExtension(input_name);
    std::string const ifile_ext = BaseLib::getFileExtension(input_name);

    std::string const obase_name = BaseLib::dropFileExtension(output_name);
    std::string const ofile_ext = BaseLib::getFileExtension(output_name);

    // Rainevents using point cloud filter
    std::cout << "Starting export... ";
    double global_min = std::numeric_limits<double>::max();
    double global_max = 0;

    std::size_t n_rasters = 0;
    for (std::size_t i = 0; ; ++i)
    {
        std::string const file_name = ibase_name + std::to_string(i) + ifile_ext;
        std::unique_ptr<GeoLib::Raster> r(
            FileIO::AsciiRasterInterface::getRasterFromASCFile(file_name));
        if (r == nullptr)
        {
            n_rasters = i;
            break;
        }

        for (auto it = r->begin(); it != r->end(); ++it)
        {
            if (*it != r->getHeader().no_data)
            {
                global_min = std::min(*it, global_min);
                global_max = std::max(*it, global_max);
            }
        }
    }

    for (std::size_t i = 0; i < n_rasters; ++i)
    {
        std::string const file_name = 
            ibase_name + std::to_string(i) + ifile_ext;
        vtkSmartPointer<VtkGeoImageSource> geo_image =
            vtkSmartPointer<VtkGeoImageSource>::New();
        geo_image->readImage(QString::fromStdString(file_name));
        VtkImageDataToPointCloudFilter *const point_cloud_filter =
            VtkImageDataToPointCloudFilter::New();
        point_cloud_filter->SetInputConnection(geo_image->GetOutputPort());
        geo_image->Update();
        point_cloud_filter->SetMinValueRange(global_min);
        point_cloud_filter->SetMaxValueRange(global_max);
        point_cloud_filter->SetMinHeight(0);
        point_cloud_filter->SetMaxHeight(max_height);
        point_cloud_filter->SetIsLinear(false);
        point_cloud_filter->SetMinNumberOfPointsPerCell(0);
        point_cloud_filter->SetMaxNumberOfPointsPerCell(max_points);
        point_cloud_filter->SetGamma(gamma);
        point_cloud_filter->Update();
        vtkSmartPointer<vtkPolyDataAlgorithm> const pd_alg =
            dynamic_cast<vtkPolyDataAlgorithm*>(point_cloud_filter);
        if (pd_alg)
        {
            std::string const output =
                obase_name + std::to_string(i) + ofile_ext;
            vtkSmartPointer<vtkXMLPolyDataWriter> pd_writer =
                vtkSmartPointer<vtkXMLPolyDataWriter>::New();
            pd_writer->SetInputData(pd_alg->GetOutputDataObject(0));
            pd_writer->SetFileName(output.c_str());
            pd_writer->Write();
        }
    }
    std::cout << "done." << std::endl;
    return EXIT_SUCCESS;
}
