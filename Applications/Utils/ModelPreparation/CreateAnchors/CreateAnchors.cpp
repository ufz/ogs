/**
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <spdlog/spdlog.h>
#include <tclap/CmdLine.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <fstream>
#include <nlohmann/json.hpp>

#include "BaseLib/FileTools.h"
#include "BaseLib/Logging.h"
#include "ComputeNaturalCoordsAlgorithm.h"
#include "InfoLib/GitInfo.h"

vtkSmartPointer<vtkUnstructuredGrid> readGrid(std::string const& input_filename)
{
    if (!BaseLib::IsFileExisting(input_filename))
    {
        OGS_FATAL("'{:s}' file does not exist.", input_filename);
    }

    vtkNew<vtkXMLUnstructuredGridReader> reader;
    reader->SetFileName(input_filename.c_str());
    reader->Update();

    if (!reader->GetOutput())
    {
        OGS_FATAL("Could not open file '{}'", input_filename);
    }
    if (reader->GetOutput()->GetNumberOfCells() == 0)
    {
        OGS_FATAL("Mesh file '{}' contains no cells.", input_filename);
    }
    return reader->GetOutput();
}

void writeGrid(vtkUnstructuredGrid* grid, std::string const& output_filename)
{
    vtkNew<vtkXMLUnstructuredGridWriter> writer;
    writer->SetFileName(output_filename.c_str());
    writer->SetInputData(grid);
    writer->Write();
}

Eigen::MatrixX3d readJSON(TCLAP::ValueArg<std::string> const& input_filename)
{
    using json = nlohmann::json;
    std::ifstream f(input_filename.getValue());
    if (!f)
    {
        OGS_FATAL("Could not open file '{:s}'", input_filename.getValue());
    }
    json const data = json::parse(f);

    auto const number_of_anchors = data["anchor_start_points"].size();
    if (number_of_anchors == 0)
    {
        OGS_FATAL("No anchors found in the json.");
    }
    if (number_of_anchors != data["anchor_end_points"].size())
    {
        OGS_FATAL(
            "Number of anchor start points does not match the number of end "
            "points");
    }

    Eigen::MatrixX3d realcoords(2 * number_of_anchors, 3);

    auto get_coordinates = [&data](const char* key, std::size_t const i)
    {
        const auto& v = data[key][i].get_ref<json::array_t const&>();
        if (v.size() != 3)
        {
            OGS_FATAL(
                "Expected a vector of length three for {:s} {}. Got vector of "
                "size {}.",
                key, i, v.size());
        }
        return Eigen::RowVector3d{v[0].get<double>(), v[1].get<double>(),
                                  v[2].get<double>()};
    };

    for (std::size_t i = 0; i < number_of_anchors; i++)
    {
        realcoords.row(2 * i).noalias() =
            get_coordinates("anchor_start_points", i);
        realcoords.row(2 * i + 1).noalias() =
            get_coordinates("anchor_end_points", i);
    }

    return realcoords;
}

int main(int argc, char** argv)
{
    // parse cmdline -----------------------------------------------------------

    TCLAP::CmdLine cmd(
        "Computes natual coordinates from given real coordinates\n\n"

        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2025, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);

    TCLAP::ValueArg<std::string> log_level_arg(
        "l", "log-level",
        "the verbosity of logging messages: none, error, warn, info, debug, "
        "all",
        false,
#ifdef NDEBUG
        "info",
#else
        "all",
#endif
        "LOG_LEVEL", cmd);

    TCLAP::ValueArg<std::string> input_filename_arg(
        "i", "input", "Input bulk mesh", true, "", "VTU_FILE", cmd);

    TCLAP::ValueArg<std::string> output_filename_arg(
        "o", "output", "Anchor mesh", true, "", "VTU_FILE", cmd);

    TCLAP::ValueArg<std::string> json_filename_arg(
        "f", "json-file", "JSON file containing anchor start and end points",
        true, "", "JSON_FILE", cmd);

    TCLAP::ValueArg<double> tolerance_arg(
        "", "tolerance", "Tolerance/Search length", false,
        std::numeric_limits<double>::epsilon(), "float", cmd);

    TCLAP::ValueArg<unsigned> max_iter_arg(
        "", "max-iter",
        "maximum number of iterations of the internal root-finding algorithm",
        false, 5, "int", cmd);

    cmd.parse(argc, argv);

    BaseLib::setConsoleLogLevel(log_level_arg.getValue());
    spdlog::set_pattern("%^%l:%$ %v");
    spdlog::set_error_handler(
        [](const std::string& msg)
        {
            std::cerr << "spdlog error: " << msg << std::endl;
            OGS_FATAL("spdlog logger error occurred.");
        });

    auto const realcoords = readJSON(json_filename_arg);
    auto const bulk_mesh = readGrid(input_filename_arg.getValue());

    // Compute natural coordinates
    namespace AU = ApplicationUtils;
    AU::ComputeNaturalCoordsResult const result = AU::computeNaturalCoords(
        bulk_mesh, realcoords, tolerance_arg.getValue(),
        max_iter_arg.getValue());
    // Write output
    auto const output_mesh = AU::toVTKGrid(result);
    writeGrid(output_mesh, output_filename_arg.getValue());

    if (!result.success)
    {
        OGS_FATAL("Failed to compute natural coordinates.");
    }

    return EXIT_SUCCESS;
}
