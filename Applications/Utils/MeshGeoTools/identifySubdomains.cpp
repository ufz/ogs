/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <tclap/CmdLine.h>

#ifdef USE_PETSC
#include <mpi.h>
#endif

#include "BaseLib/RunTime.h"
#include "InfoLib/GitInfo.h"
#include "MeshGeoToolsLib/IdentifySubdomainMesh.h"
#include "MeshGeoToolsLib/MeshNodeSearcher.h"
#include "MeshGeoToolsLib/SearchLength.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/IO/writeMeshToFile.h"
#include "MeshLib/Mesh.h"

std::vector<std::unique_ptr<MeshLib::Mesh>> readMeshes(
    std::vector<std::string> const& filenames)
{
    std::vector<std::unique_ptr<MeshLib::Mesh>> meshes;
    meshes.reserve(filenames.size());

    for (auto const& filename : filenames)
    {
        auto mesh = MeshLib::IO::readMeshFromFile(filename);
        if (mesh == nullptr)
        {
            OGS_FATAL("Could not read mesh from '{:s}' file.", filename);
        }
        meshes.emplace_back(mesh);
    }
    if (meshes.empty())
    {
        OGS_FATAL("No subdomain meshes were read.");
    }
    return meshes;
}

int main(int argc, char* argv[])
{
    BaseLib::RunTime run_time;
    run_time.start();

    TCLAP::CmdLine cmd(
        "Checks if the subdomain meshes are part of the bulk mesh and writes "
        "the 'bulk_node_ids' and the 'bulk_element_ids' in each of them. The "
        "documentation is available at "
        "https://www.opengeosys.org/docs/tools/meshing-submeshes/"
        "identifysubdomains/.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2024, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);

    TCLAP::SwitchArg force_overwrite_arg(
        "f", "force", "Overwriting existing subdomain meshes.");
    cmd.add(force_overwrite_arg);

    TCLAP::ValueArg<std::string> output_prefix_arg(
        "o",
        "output_prefix",
        "Prefix the subdomain meshes' filenames with the output prefix/path.",
        false,
        "",
        "path");
    cmd.add(output_prefix_arg);

    TCLAP::ValueArg<double> search_length_arg(
        "s",
        "searchlength",
        "search length determining radius for the node search algorithm. "
        "Non-negative floating point number (default 1e-16) ",
        false,
        1e-16,
        "float");
    cmd.add(search_length_arg);

    TCLAP::ValueArg<std::string> bulk_mesh_arg(
        "m", "mesh", "the file name of the bulk mesh", true, "", "mesh file");
    cmd.add(bulk_mesh_arg);

    // All the remaining arguments are used as file names for boundary/subdomain
    // meshes.
    TCLAP::UnlabeledMultiArg<std::string> subdomain_meshes_filenames_arg(
        "subdomain_meshes_filenames", "mesh file names.", true,
        "subdomain mesh file");
    cmd.add(subdomain_meshes_filenames_arg);
    cmd.parse(argc, argv);

#ifdef USE_PETSC
    MPI_Init(&argc, &argv);
#endif

    //
    // The bulk mesh.
    //
    BaseLib::RunTime mesh_reading_time;
    mesh_reading_time.start();
    std::unique_ptr<MeshLib::Mesh> bulk_mesh{
        MeshLib::IO::readMeshFromFile(bulk_mesh_arg.getValue())};
    if (bulk_mesh == nullptr)
    {
        OGS_FATAL("Could not read bulk mesh from '{:s}'",
                  bulk_mesh_arg.getValue());
    }

    //
    // Read the subdomain meshes.
    //
    auto const subdomain_meshes =
        readMeshes(subdomain_meshes_filenames_arg.getValue());
    INFO("Mesh reading time: {:g} s", mesh_reading_time.elapsed());

    //
    // Bulk mesh node searcher.
    //
    BaseLib::RunTime mesh_node_searcher_construction_time;
    mesh_node_searcher_construction_time.start();
    auto const& mesh_node_searcher =
        MeshGeoToolsLib::MeshNodeSearcher::getMeshNodeSearcher(
            *bulk_mesh,
            std::make_unique<MeshGeoToolsLib::SearchLength>(
                search_length_arg.getValue()));
    INFO("MeshNodeSearcher construction time: {:g} s",
         mesh_node_searcher_construction_time.elapsed());

    //
    // Identify the subdomains in the bulk mesh.
    //
    BaseLib::RunTime identify_subdomain_time;
    identify_subdomain_time.start();
    for (auto& mesh_ptr : subdomain_meshes)
    {
        // If force overwrite is set or the output is to different mesh than
        // the input mesh.
        bool const overwrite_property_vectors =
            force_overwrite_arg.getValue() ||
            !output_prefix_arg.getValue().empty();
        identifySubdomainMesh(*mesh_ptr, *bulk_mesh, mesh_node_searcher,
                              overwrite_property_vectors);
    }
    INFO("identifySubdomains time: {:g} s", identify_subdomain_time.elapsed());

    //
    // Output after the successful subdomain mesh identification.
    //
    BaseLib::RunTime writing_time;
    writing_time.start();
    for (auto const& mesh_ptr : subdomain_meshes)
    {
        MeshLib::IO::writeMeshToFile(
            *mesh_ptr,
            output_prefix_arg.getValue() + mesh_ptr->getName() + ".vtu");
    }
    INFO("writing time: {:g} s", writing_time.elapsed());

#ifdef USE_PETSC
    MPI_Finalize();
#endif
    INFO("Entire run time: {:g} s", run_time.elapsed());
    return EXIT_SUCCESS;
}
