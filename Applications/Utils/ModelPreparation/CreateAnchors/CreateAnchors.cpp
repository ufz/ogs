/**
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <tclap/CmdLine.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <fstream>
#include <nlohmann/json.hpp>

#include "BaseLib/FileTools.h"
#include "BaseLib/Logging.h"
#include "BaseLib/MPI.h"
#include "BaseLib/TCLAPArguments.h"
#include "ComputeNaturalCoordsAlgorithm.h"
#include "InfoLib/GitInfo.h"
#include "MeshLib/IO/VtkIO/VtuInterface.h"

namespace AU = ApplicationUtils;

vtkSmartPointer<vtkUnstructuredGrid> readGrid(std::string const& input_filename)
{
    vtkSmartPointer<vtkUnstructuredGrid> grid =
        MeshLib::IO::VtuInterface::readVtuFileToVtkUnstructuredGrid(
            input_filename);
    if (grid == nullptr)
    {
        OGS_FATAL("Could not open file '{}'", input_filename);
    }
    if (grid->GetNumberOfCells() == 0)
    {
        OGS_FATAL("Mesh file '{}' contains no cells.", input_filename);
    }
    return grid;
}

void writeGrid(vtkUnstructuredGrid* grid, std::string const& output_filename)
{
    vtkNew<vtkXMLUnstructuredGridWriter> writer;
    writer->SetFileName(output_filename.c_str());
    writer->SetInputData(grid);
    writer->Write();
}

void checkJSONEntries(nlohmann::json const& data, size_t number_of_anchors)
{
    std::vector<std::string> required_keys = {"anchor_stiffness",
                                              "anchor_cross_sectional_area"};
    std::vector<std::string> optional_keys = {"maximum_anchor_stress",
                                              "initial_anchor_stress",
                                              "residual_anchor_stress"};
    if (number_of_anchors == 0)
    {
        OGS_FATAL("No anchors found in the json.");
    }
    for (const auto& key : required_keys)
    {
        if (!data.contains(key))
        {
            OGS_FATAL("JSON file does not contain required key '{:s}'", key);
        }
        if (number_of_anchors != data[key].size())
        {
            OGS_FATAL(
                "Number of anchor start points does not match the number of {} "
                "entries.",
                key);
        }
        for (size_t i = 0; i < number_of_anchors; ++i)
        {
            if (!data[key][i].is_number())
            {
                OGS_FATAL("Non-numeric element in JSON array for key {}!", key);
            }
        }
    }
    for (const auto& key : optional_keys)
    {
        if (data.contains(key))
        {
            if (number_of_anchors != data[key].size())
            {
                OGS_FATAL(
                    "Number of anchor start points does not match the number "
                    "of {} "
                    "entries.",
                    key);
            }
            for (size_t i = 0; i < number_of_anchors; ++i)
            {
                if (!data[key][i].is_number())
                {
                    OGS_FATAL("Non-numeric element in JSON array for key {}!",
                              key);
                }
            }
        }
    }
}

AU::ComputeNaturalCoordsResult readJSON(
    TCLAP::ValueArg<std::string> const& input_filename)
{
    using json = nlohmann::json;
    std::ifstream f(input_filename.getValue());
    if (!f)
    {
        OGS_FATAL("Could not open file '{:s}'", input_filename.getValue());
    }
    json const data = json::parse(f);

    auto const number_of_anchors = data["anchor_start_points"].size();
    if (number_of_anchors != data["anchor_end_points"].size())
    {
        OGS_FATAL(
            "Number of anchor start points does not match the number of "
            "anchor end points.");
    }
    checkJSONEntries(data, number_of_anchors);

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
    Eigen::VectorXd initial_anchor_stress(number_of_anchors);
    Eigen::VectorXd maximum_anchor_stress(number_of_anchors);
    Eigen::VectorXd residual_anchor_stress(number_of_anchors);
    Eigen::VectorXd anchor_cross_sectional_area(number_of_anchors);
    Eigen::VectorXd anchor_stiffness(number_of_anchors);
    if (!data.contains("initial_anchor_stress"))
    {
        initial_anchor_stress.setZero();
    }
    else
    {
        for (size_t i = 0; i < number_of_anchors; ++i)
        {
            initial_anchor_stress(i) =
                data["initial_anchor_stress"][i].get<double>();
        }
    }
    if (!data.contains("maximum_anchor_stress"))
    {
        maximum_anchor_stress = Eigen::VectorXd::Constant(
            number_of_anchors, std::numeric_limits<double>::max());
    }
    else
    {
        for (size_t i = 0; i < number_of_anchors; ++i)
        {
            maximum_anchor_stress(i) =
                data["maximum_anchor_stress"][i].get<double>();
        }
    }
    if (!data.contains("residual_anchor_stress"))
    {
        residual_anchor_stress = Eigen::VectorXd::Constant(
            number_of_anchors, std::numeric_limits<double>::max());
    }
    else
    {
        for (size_t i = 0; i < number_of_anchors; ++i)
        {
            residual_anchor_stress(i) =
                data["residual_anchor_stress"][i].get<double>();
        }
    }
    for (size_t i = 0; i < number_of_anchors; ++i)
    {
        anchor_cross_sectional_area(i) =
            data["anchor_cross_sectional_area"][i].get<double>();
        anchor_stiffness(i) = data["anchor_stiffness"][i].get<double>();
    }
    AU::ComputeNaturalCoordsResult json_data;
    json_data.real_coords = realcoords;
    json_data.initial_anchor_stress = initial_anchor_stress;
    json_data.maximum_anchor_stress = maximum_anchor_stress;
    json_data.residual_anchor_stress = residual_anchor_stress;
    json_data.anchor_cross_sectional_area = anchor_cross_sectional_area;
    json_data.anchor_stiffness = anchor_stiffness;
    return json_data;
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

    auto log_level_arg = BaseLib::makeLogLevelArg();
    cmd.add(log_level_arg);

    TCLAP::ValueArg<std::string> input_filename_arg(
        "i", "input", "Input (.vtu) bulk mesh file", true, "", "INPUT_FILE",
        cmd);

    TCLAP::ValueArg<std::string> output_filename_arg(
        "o", "output", "Output (.vtu). Anchor mesh file", true, "",
        "OUTPUT_FILE", cmd);

    TCLAP::ValueArg<std::string> json_filename_arg(
        "f", "json-file",
        "Input (.json) JSON file containing anchor start and end points", true,
        "", "INPUT_FILE", cmd);

    TCLAP::ValueArg<double> tolerance_arg(
        "", "tolerance", "Tolerance/Search length", false,
        std::numeric_limits<double>::epsilon(), "float", cmd);

    TCLAP::ValueArg<unsigned> max_iter_arg(
        "", "max-iter",
        "maximum number of iterations of the internal root-finding "
        "algorithm, (min = 0)",
        false, 5, "MAX_ITER", cmd);

    cmd.parse(argc, argv);

    BaseLib::MPI::Setup mpi_setup(argc, argv);
    BaseLib::initOGSLogger(log_level_arg.getValue());

    AU::ComputeNaturalCoordsResult const json_data =
        readJSON(json_filename_arg);
    auto const bulk_mesh = readGrid(input_filename_arg.getValue());

    // Compute natural coordinates
    AU::ComputeNaturalCoordsResult const result =
        AU::computeNaturalCoords(bulk_mesh, json_data, tolerance_arg.getValue(),
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
