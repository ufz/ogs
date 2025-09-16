/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <tclap/CmdLine.h>

#include <fstream>
#include <range/v3/algorithm/minmax.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/filter.hpp>
#include <range/v3/view/iota.hpp>

#include "BaseLib/FileTools.h"
#include "BaseLib/Logging.h"
#include "BaseLib/MPI.h"
#include "BaseLib/TCLAPArguments.h"
#include "InfoLib/GitInfo.h"
#include "MeshLib/IO/VtkIO/VtuInterface.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/Mesh.h"
#include "MeshToolsLib/MeshEditing/RemoveMeshComponents.h"

MeshLib::Mesh* extractMatGroup(MeshLib::Mesh const& mesh, int const mat_id)
{
    auto const mat_ids =
        *mesh.getProperties().getPropertyVector<int>("MaterialIDs");

    auto const elem_list =
        ranges::views::iota(std::size_t{0}, mat_ids.size()) |
        ranges::views::filter([&](std::size_t i)
                              { return mat_ids[i] != mat_id; }) |
        ranges::to<std::vector>;

    return MeshToolsLib::removeElements(mesh, elem_list, "matgroup");
}

int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Takes a mesh with multiple MaterialIDs and writes elements of a given "
        "ID into a new mesh. If no ID is specified, meshes for each existing "
        "MaterialID are created.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2025, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);
    TCLAP::ValueArg<std::size_t> arg_mat_id(
        "m", "material-id",
        "The MaterialID for which elements should be extracted into a new "
        "mesh.",
        false, 0, "MATERIAL_ID");
    cmd.add(arg_mat_id);
    TCLAP::ValueArg<std::string> output_arg(
        "o", "output", "Output (.vtu). Name of the output mesh", true, "",
        "OUTPUT_FILE");
    cmd.add(output_arg);
    TCLAP::ValueArg<std::string> input_arg(
        "i", "input", "Input (.vtu). Name of the input mesh", true, "",
        "INPUT_FILE");
    cmd.add(input_arg);
    auto log_level_arg = BaseLib::makeLogLevelArg();
    cmd.add(log_level_arg);
    cmd.parse(argc, argv);

    BaseLib::MPI::Setup mpi_setup(argc, argv);
    BaseLib::initOGSLogger(log_level_arg.getValue());

    std::string const input_name = input_arg.getValue();
    std::string const output_name = output_arg.getValue();
    std::string const base_name = BaseLib::dropFileExtension(output_name);
    std::string const ext = BaseLib::getFileExtension(output_name);

    std::unique_ptr<MeshLib::Mesh> const mesh(
        MeshLib::IO::readMeshFromFile(input_name));
    if (mesh == nullptr)
    {
        ERR("Error reading input mesh. Aborting...");
        return EXIT_FAILURE;
    }

    auto mat_ids = MeshLib::materialIDs(*mesh);
    if (mat_ids == nullptr)
    {
        ERR("No material IDs found in mesh. Aborting...");
        return EXIT_FAILURE;
    }

    auto const [min, max] = ranges::minmax(*mat_ids);
    if (min == max)
    {
        ERR("Mesh only contains one material, no extraction required.");
        return EXIT_FAILURE;
    }
    int min_id, max_id;
    if (arg_mat_id.isSet())
    {
        min_id = static_cast<int>(arg_mat_id.getValue());
        if (min_id < min || min_id > max)
        {
            ERR("Specified material ID does not exist.");
            return EXIT_FAILURE;
        }
        max_id = min_id;
    }
    else
    {
        min_id = min;
        max_id = max;
    }

    std::ofstream ostream;
    if (min_id != max_id)
    {
        ostream.open(base_name + "_Layers.txt");
    }

    for (int i = min_id; i <= max_id; ++i)
    {
        INFO("Extracting material group {:d}...", i);
        std::unique_ptr<MeshLib::Mesh> mat_group(extractMatGroup(*mesh, i));
        if (mat_group == nullptr)
        {
            WARN("No elements with material group {:d} found.", i);
            continue;
        }
        MeshLib::IO::VtuInterface vtu(mat_group.get());
        std::string const file_name(base_name + "_Layer" + std::to_string(i) +
                                    ext);
        vtu.writeToFile(file_name);
        if (ostream.is_open())
        {
            ostream << file_name << "\n";
        }
    }
    if (ostream.is_open())
    {
        ostream.close();
    }
    return EXIT_SUCCESS;
}
