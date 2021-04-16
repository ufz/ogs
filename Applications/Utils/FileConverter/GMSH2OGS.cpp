/**
 * \file
 * \author Thomas Fischer
 * \date   2011-12-13
 * \brief  Implementation of the GMSH2OGS converter.
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// STL
#include <algorithm>
#include <string>

// ThirdParty
#include <tclap/CmdLine.h>

// BaseLib
#include "BaseLib/FileTools.h"
#include "BaseLib/RunTime.h"
#include "InfoLib/GitInfo.h"
#ifndef WIN32
#include "BaseLib/MemWatch.h"
#endif

#include "Applications/FileIO/Gmsh/GmshReader.h"
#include "GeoLib/AABB.h"
#include "MeshGeoToolsLib/IdentifySubdomainMesh.h"
#include "MeshGeoToolsLib/MeshNodeSearcher.h"
#include "MeshLib/IO/writeMeshToFile.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshEditing/DuplicateMeshComponents.h"
#include "MeshLib/MeshEditing/RemoveMeshComponents.h"
#include "MeshLib/MeshInformation.h"
#include "MeshLib/MeshQuality/MeshValidation.h"
#include "MeshLib/MeshSearch/ElementSearch.h"

static std::unique_ptr<MeshLib::Mesh> createMeshFromElements(
    MeshLib::Mesh const& mesh,
    std::vector<MeshLib::Element*> const& selected_elements,
    std::string mesh_name)
{
    // Create boundary mesh.
    auto nodes = copyNodeVector(mesh.getNodes());
    auto elements = copyElementVector(selected_elements, nodes);

    // Cleanup unused nodes
    removeMarkedNodes(markUnusedNodes(elements, nodes), nodes);

    return std::make_unique<MeshLib::Mesh>(std::move(mesh_name), nodes,
                                           elements);
}

static std::vector<std::unique_ptr<MeshLib::Mesh>> extractBoundaryMeshes(
    MeshLib::Mesh const& mesh, std::vector<std::size_t> selected_element_ids)
{
    auto const material_ids = materialIDs(mesh);
    if (material_ids == nullptr)
    {
        OGS_FATAL(
            "GMSH2OGS: Expected material ids to be present in the mesh to "
            "extract boundary meshes.");
    }

    std::vector<std::unique_ptr<MeshLib::Mesh>> boundary_meshes;

    auto const& elements = mesh.getElements();

    while (!selected_element_ids.empty())
    {
        // Partition in two blocks, with elements for the material id at
        // the end, s.t. one can erase them easily.
        int const material_id = (*material_ids)[selected_element_ids.back()];
        auto split = std::partition(
            begin(selected_element_ids), end(selected_element_ids),
            [&material_id, &material_ids](int const id) {
                return (*material_ids)[id] != material_id;
            });

        // Add elements with same material id to the mesh.
        std::vector<MeshLib::Element*> single_material_elements;
        single_material_elements.reserve(
            std::distance(split, end(selected_element_ids)));
        std::transform(split, end(selected_element_ids),
                       std::back_inserter(single_material_elements),
                       [&](int const id) { return elements[id]; });

        // Remove already extracted elements.
        selected_element_ids.erase(split, end(selected_element_ids));

        // Create boundary mesh and identify the nodes/elements.
        boundary_meshes.push_back(createMeshFromElements(
            mesh, single_material_elements, std::to_string(material_id)));
    }
    return boundary_meshes;
}

static void identifyAndWriteBoundaryMeshes(
    MeshLib::Mesh const& mesh,
    std::string const& file_name,
    std::vector<std::unique_ptr<MeshLib::Mesh>>& boundary_meshes)
{
    // Bulk mesh node searcher usef for boundary mesh identification.
    auto const& mesh_node_searcher =
        MeshGeoToolsLib::MeshNodeSearcher::getMeshNodeSearcher(
            mesh,
            std::make_unique<MeshGeoToolsLib::SearchLength>(
                0));  // Exact match of nodes.

    for (auto& boundary_mesh : boundary_meshes)
    {
        identifySubdomainMesh(*boundary_mesh, mesh, mesh_node_searcher);

        // Save the boundary mesh.
        auto boundary_mesh_file_name = BaseLib::dropFileExtension(file_name) +
                                       '_' + boundary_mesh->getName() +
                                       BaseLib::getFileExtension(file_name);

        MeshLib::IO::writeMeshToFile(*boundary_mesh, boundary_mesh_file_name);
    }
}

int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Converting meshes in gmsh file format (ASCII, version 2.2) to a vtk "
        "unstructured grid file (new OGS file format) or to the old OGS file "
        "format - see options.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2021, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);

    TCLAP::ValueArg<std::string> ogs_mesh_arg(
        "o",
        "out",
        "filename for output mesh (if extension is .msh, old OGS-5 fileformat "
        "is written, if extension is .vtu, a vtk unstructure grid file is "
        "written (OGS-6 mesh format))",
        true,
        "",
        "filename as string");
    cmd.add(ogs_mesh_arg);

    TCLAP::ValueArg<std::string> gmsh_mesh_arg("i", "in", "gmsh input file",
                                               true, "", "filename as string");
    cmd.add(gmsh_mesh_arg);

    TCLAP::SwitchArg valid_arg("v", "validation", "validate the mesh");
    cmd.add(valid_arg);

    TCLAP::SwitchArg create_boundary_meshes_arg(
        "b", "boundaries", "if set, boundary meshes will be generated");
    cmd.add(create_boundary_meshes_arg);

    TCLAP::SwitchArg exclude_lines_arg(
        "e", "exclude-lines",
        "if set, lines will not be written to the ogs mesh");
    cmd.add(exclude_lines_arg);

    cmd.parse(argc, argv);

    // *** read mesh
    INFO("Reading {:s}.", gmsh_mesh_arg.getValue());
#ifndef WIN32
    BaseLib::MemWatch mem_watch;
    unsigned long mem_without_mesh(mem_watch.getVirtMemUsage());
#endif
    BaseLib::RunTime run_time;
    run_time.start();
    MeshLib::Mesh* mesh(FileIO::GMSH::readGMSHMesh(gmsh_mesh_arg.getValue()));

    if (mesh == nullptr)
    {
        INFO("Could not read mesh from {:s}.", gmsh_mesh_arg.getValue());
        return -1;
    }
#ifndef WIN32
    INFO("Mem for mesh: {} MiB",
         (mem_watch.getVirtMemUsage() - mem_without_mesh) / (1024 * 1024));
#endif

    INFO("Time for reading: {:f} seconds.", run_time.elapsed());
    INFO("Read {:d} nodes and {:d} elements.", mesh->getNumberOfNodes(),
         mesh->getNumberOfElements());

    // Optionally remove line elements or create boundary meshes.
    if (exclude_lines_arg.getValue() || create_boundary_meshes_arg.getValue())
    {
        auto ex = MeshLib::ElementSearch(*mesh);
        ex.searchByElementType(MeshLib::MeshElemType::LINE);
        auto const& selected_element_ids = ex.getSearchedElementIDs();

        // First we extract the boundary meshes, then optionally remove the line
        // elements, and only then run the node/element identification and write
        // the meshes.

        std::vector<std::unique_ptr<MeshLib::Mesh>> boundary_meshes;
        if (create_boundary_meshes_arg.getValue())
        {
            boundary_meshes =
                extractBoundaryMeshes(*mesh, selected_element_ids);
        }

        if (exclude_lines_arg.getValue())
        {
            auto m = MeshLib::removeElements(*mesh, selected_element_ids,
                                             mesh->getName() + "-withoutLines");
            if (m != nullptr)
            {
                INFO("Removed {:d} lines.",
                     mesh->getNumberOfElements() - m->getNumberOfElements());
                std::swap(m, mesh);
                delete m;
            }
            else
            {
                INFO("Mesh does not contain any lines.");
            }
        }

        if (create_boundary_meshes_arg.getValue())
        {
            identifyAndWriteBoundaryMeshes(*mesh, ogs_mesh_arg.getValue(),
                                           boundary_meshes);
        }
    }
    // *** print meshinformation

    INFO("Please check your mesh carefully!");
    INFO(
        "Degenerated or redundant mesh elements can cause OGS to stop or "
        "misbehave.");
    INFO("Use the -e option to delete redundant line elements.");

    // Geometric information
    const GeoLib::AABB aabb = MeshLib::MeshInformation::getBoundingBox(*mesh);
    auto const minPt(aabb.getMinPoint());
    auto const maxPt(aabb.getMaxPoint());
    INFO("Node coordinates:");
    INFO("\tx [{:g}, {:g}] (extent {:g})", minPt[0], maxPt[0],
         maxPt[0] - minPt[0]);
    INFO("\ty [{:g}, {:g}] (extent {:g})", minPt[1], maxPt[1],
         maxPt[1] - minPt[1]);
    INFO("\tz [{:g}, {:g}] (extent {:g})", minPt[2], maxPt[2],
         maxPt[2] - minPt[2]);

    INFO("Edge length: [{:g}, {:g}]", mesh->getMinEdgeLength(),
         mesh->getMaxEdgeLength());

    // Element information
    MeshLib::MeshInformation::writeAllNumbersOfElementTypes(*mesh);

    MeshLib::MeshInformation::writePropertyVectorInformation(*mesh);

    if (valid_arg.isSet())
    {
        MeshLib::MeshInformation::writeMeshValidationResults(*mesh);
    }

    // *** write mesh in new format
    MeshLib::IO::writeMeshToFile(*mesh, ogs_mesh_arg.getValue());

    delete mesh;
}
