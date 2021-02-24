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
#include <string>
#include <algorithm>

// ThirdParty
#include <tclap/CmdLine.h>

// BaseLib
#include "InfoLib/GitInfo.h"
#include "BaseLib/FileTools.h"
#include "BaseLib/RunTime.h"
#ifndef WIN32
#include "BaseLib/MemWatch.h"
#endif

#include "Applications/FileIO/Gmsh/GmshReader.h"
#include "GeoLib/AABB.h"
#include "MeshLib/IO/writeMeshToFile.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshEditing/DuplicateMeshComponents.h"
#include "MeshLib/MeshEditing/RemoveMeshComponents.h"
#include "MeshLib/MeshInformation.h"
#include "MeshLib/MeshQuality/MeshValidation.h"
#include "MeshLib/MeshSearch/ElementSearch.h"

static MeshLib::Mesh createMeshFromElements(
    MeshLib::Mesh const& mesh,
    std::vector<MeshLib::Element*> const& selected_elements)
{
    // Create boundary mesh.
    auto nodes = copyNodeVector(mesh.getNodes());
    auto elements = copyElementVector(selected_elements, nodes);

    // Cleanup unused nodes
    removeMarkedNodes(markUnusedNodes(elements, nodes), nodes);

    return {mesh.getName(), nodes, elements};
}

static void extractBoundaryMeshes(MeshLib::Mesh const& mesh,
                                  MeshLib::ElementSearch const& ex,
                                  std::string const& file_name)
{
    // Copy line element ids to process them to boundary meshes.
    auto line_element_ids = ex.getSearchedElementIDs();

    auto const material_ids = materialIDs(mesh);
    if (material_ids == nullptr)
    {
        OGS_FATAL(
            "GMSH2OGS: Expected material ids to be present in the mesh to "
            "extract boundary meshes.");
    }
    auto const& elements = mesh.getElements();

    while (!line_element_ids.empty())
    {
        // Partition in two blocks, with elements for the material id at
        // the end, s.t. one can erase them easily.
        int const material_id = (*material_ids)[line_element_ids.back()];
        auto split =
            std::partition(begin(line_element_ids), end(line_element_ids),
                           [&material_id, &material_ids](int const id) {
                               return (*material_ids)[id] != material_id;
                           });

        // Add elements with same material id to the mesh.
        std::vector<MeshLib::Element*> single_material_elements;
        single_material_elements.reserve(
            std::distance(split, end(line_element_ids)));
        std::transform(split, end(line_element_ids),
                       std::back_inserter(single_material_elements),
                       [&](int const id) { return elements[id]; });

        // Remove already extracted elements.
        line_element_ids.erase(split, end(line_element_ids));

        // Create boundary mesh and save.
        auto boundary_mesh =
            createMeshFromElements(mesh, single_material_elements);

        auto boundary_mesh_file_name = BaseLib::dropFileExtension(file_name) +
                                       '_' + std::to_string(material_id) +
                                       BaseLib::getFileExtension(file_name);

        MeshLib::IO::writeMeshToFile(boundary_mesh, boundary_mesh_file_name);
    }
}

int main (int argc, char* argv[])
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
        "filename for output mesh (if extension is .msh, old OGS-5 fileformat is written, if extension is .vtu, a vtk unstructure grid file is written (OGS-6 mesh format))",
        true,
        "",
        "filename as string");
    cmd.add(ogs_mesh_arg);

    TCLAP::ValueArg<std::string> gmsh_mesh_arg(
        "i",
        "in",
        "gmsh input file",
        true,
        "",
        "filename as string");
    cmd.add(gmsh_mesh_arg);

    TCLAP::SwitchArg valid_arg("v", "validation", "validate the mesh");
    cmd.add(valid_arg);

    TCLAP::SwitchArg create_boundary_meshes_arg(
        "b", "boundaries", "if set, boundary meshes will be generated");
    cmd.add(create_boundary_meshes_arg);

    TCLAP::SwitchArg exclude_lines_arg("e", "exclude-lines",
        "if set, lines will not be written to the ogs mesh");
    cmd.add(exclude_lines_arg);

    cmd.parse(argc, argv);

    // *** read mesh
    INFO("Reading {:s}.", gmsh_mesh_arg.getValue());
#ifndef WIN32
    BaseLib::MemWatch mem_watch;
    unsigned long mem_without_mesh (mem_watch.getVirtMemUsage());
#endif
    BaseLib::RunTime run_time;
    run_time.start();
    MeshLib::Mesh* mesh(
        FileIO::GMSH::readGMSHMesh(gmsh_mesh_arg.getValue()));

    if (mesh == nullptr) {
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

    // *** remove line elements on request
    if (exclude_lines_arg.getValue() || create_boundary_meshes_arg.getValue())
    {
        auto ex = MeshLib::ElementSearch(*mesh);
        ex.searchByElementType(MeshLib::MeshElemType::LINE);
        if (create_boundary_meshes_arg.getValue())
        {
            extractBoundaryMeshes(*mesh, ex, ogs_mesh_arg.getValue());
        }
        if (exclude_lines_arg.getValue())
        {
            auto m = MeshLib::removeElements(*mesh, ex.getSearchedElementIDs(),
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

