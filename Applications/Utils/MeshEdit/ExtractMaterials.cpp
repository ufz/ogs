/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <fstream>

// ThirdParty
#include <tclap/CmdLine.h>

#include "InfoLib/GitInfo.h"

#include "BaseLib/FileTools.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/IO/VtkIO/VtuInterface.h"
#include "MeshLib/MeshEditing/RemoveMeshComponents.h"

MeshLib::Mesh* extractMatGroup(MeshLib::Mesh const& mesh, int const mat_id)
{
    std::vector<std::size_t> elem_list;
    std::vector<int> const mat_ids =
        *mesh.getProperties().getPropertyVector<int>("MaterialIDs");
    std::size_t const n_elems = mat_ids.size();
    for (std::size_t i = 0; i < n_elems; ++i)
    {
        if (mat_ids[i] != mat_id)
        {
            elem_list.push_back(i);
        }
    }

    if (elem_list.empty())
    {
        return nullptr;
    }

    return MeshLib::removeElements(mesh, elem_list, "matgroup");
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
            "Copyright (c) 2012-2021, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);
    TCLAP::ValueArg<std::size_t> arg_mat_id(
        "m", "material-id",
        "The MaterialID for which elements should be extracted into a new "
        "mesh.",
        false, 0, "Number specifying the MaterialID");
    cmd.add(arg_mat_id);
    TCLAP::ValueArg<std::string> output_arg(
        "o", "output",
        "Name of the output mesh (*.vtu)",
        true, "", "output file name");
    cmd.add(output_arg);
    TCLAP::ValueArg<std::string> input_arg(
        "i", "input",
        "Name of the input mesh (*.vtu)",
        true, "", "input file name");
    cmd.add(input_arg);
    cmd.parse( argc, argv );

    std::string const input_name = input_arg.getValue();
    std::string const output_name = output_arg.getValue();
    std::string const base_name = BaseLib::dropFileExtension(output_name);
    std::string const ext = BaseLib::getFileExtension(output_name);

    std::unique_ptr<MeshLib::Mesh> const mesh (MeshLib::IO::readMeshFromFile(input_name));
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

    auto id_range = std::minmax_element(mat_ids->cbegin(), mat_ids->cend());
    if (id_range.first == id_range.second)
    {
        ERR("Mesh only contains one material, no extraction required.");
        return EXIT_FAILURE;
    }
    int min_id, max_id;
    if (arg_mat_id.isSet())
    {
        min_id = static_cast<int>(arg_mat_id.getValue());
        if (min_id < *id_range.first || min_id > *id_range.second)
        {
            ERR("Specified material ID does not exist.");
            return EXIT_FAILURE;
        }
        max_id = min_id;
    }
    else
    {
        min_id = *id_range.first;
        max_id = *id_range.second;
    }

    std::ofstream ostream;
    if (min_id != max_id)
    {
        ostream.open(base_name + "_Layers.txt");
    }

    for (int i = min_id; i <= max_id; ++i)
    {
        INFO("Extracting material group {:d}...", i);
        std::unique_ptr<MeshLib::Mesh> mat_group (extractMatGroup(*mesh, i));
        if (mat_group == nullptr)
        {
            WARN("No elements with material group {:d} found.", i);
            continue;
        }
        MeshLib::IO::VtuInterface vtu(mat_group.get());
        std::string const file_name(base_name + "_Layer" + std::to_string(i) + ext);
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
