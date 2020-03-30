/**
 * \file
 * \author Thomas Fischer
 * \date   2011-12-13
 * \brief  Implementation of the GMSH2OGS converter.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
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

#include "Applications/ApplicationsLib/LogogSetup.h"

// BaseLib
#include "InfoLib/GitInfo.h"
#include "BaseLib/FileTools.h"
#include "BaseLib/RunTime.h"
#ifndef WIN32
#include "BaseLib/MemWatch.h"
#endif

#include "GeoLib/AABB.h"

#include "Applications/FileIO/Gmsh/GmshReader.h"

#include "MeshLib/IO/writeMeshToFile.h"
#include "MeshLib/MeshSearch/ElementSearch.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshEditing/RemoveMeshComponents.h"
#include "MeshLib/MeshInformation.h"
#include "MeshLib/MeshQuality/MeshValidation.h"

int main (int argc, char* argv[])
{
    ApplicationsLib::LogogSetup logog_setup;

    TCLAP::CmdLine cmd(
        "Converting meshes in gmsh file format (ASCII, version 2.2) to a vtk "
        "unstructured grid file (new OGS file format) or to the old OGS file "
        "format - see options.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2020, OpenGeoSys Community "
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

    TCLAP::SwitchArg valid_arg("v","validation","validate the mesh");
    cmd.add( valid_arg );

    TCLAP::SwitchArg exclude_lines_arg("e", "exclude-lines",
        "if set, lines will not be written to the ogs mesh");
    cmd.add(exclude_lines_arg);

    cmd.parse(argc, argv);

    // *** read mesh
    INFO("Reading %s.", gmsh_mesh_arg.getValue().c_str());
#ifndef WIN32
    BaseLib::MemWatch mem_watch;
    unsigned long mem_without_mesh (mem_watch.getVirtMemUsage());
#endif
    BaseLib::RunTime run_time;
    run_time.start();
    MeshLib::Mesh* mesh(
        FileIO::GMSH::readGMSHMesh(gmsh_mesh_arg.getValue()));

    if (mesh == nullptr) {
        INFO("Could not read mesh from %s.", gmsh_mesh_arg.getValue().c_str());
        return -1;
    }
#ifndef WIN32
    INFO("Mem for mesh: %i MB",
         (mem_watch.getVirtMemUsage() - mem_without_mesh) / (1024 * 1024));
#endif

    INFO("Time for reading: %f seconds.", run_time.elapsed());
    INFO("Read %d nodes and %d elements.", mesh->getNumberOfNodes(), mesh->getNumberOfElements());

    // *** remove line elements on request
    if (exclude_lines_arg.getValue()) {
        auto ex = MeshLib::ElementSearch(*mesh);
        ex.searchByElementType(MeshLib::MeshElemType::LINE);
        auto m = MeshLib::removeElements(*mesh, ex.getSearchedElementIDs(), mesh->getName()+"-withoutLines");
        if (m != nullptr) {
            INFO("Removed %d lines.", mesh->getNumberOfElements() - m->getNumberOfElements());
            std::swap(m, mesh);
            delete m;
        } else {
            INFO("Mesh does not contain any lines.");
        }
    }
    // *** print meshinformation
    INFO("Please check your mesh carefully!");
    INFO("Degenerated or redundant mesh elements can cause OGS to stop or misbehave.");
    INFO("Use the -e option to delete redundant line elements.");
    // Geometric information
    const GeoLib::AABB aabb(MeshLib::MeshInformation::getBoundingBox(*mesh));
    auto minPt(aabb.getMinPoint());
    auto maxPt(aabb.getMaxPoint());
    INFO("Node coordinates:");
    INFO("\tx [%g, %g] (extent %g)", minPt[0], maxPt[0], maxPt[0]-minPt[0]);
    INFO("\ty [%g, %g] (extent %g)", minPt[1], maxPt[1], maxPt[1]-minPt[1]);
    INFO("\tz [%g, %g] (extent %g)", minPt[2], maxPt[2], maxPt[2]-minPt[2]);

    INFO("Edge length: [%g, %g]", mesh->getMinEdgeLength(), mesh->getMaxEdgeLength());

    // Element information
    const std::array<unsigned, 7> nr_ele_types(MeshLib::MeshInformation::getNumberOfElementTypes(*mesh));
    INFO("Element types:");
    unsigned etype = 0;
    if (nr_ele_types[etype] > 0)
    {
        INFO("\t%d lines", nr_ele_types[etype]);
    }
    if (nr_ele_types[++etype] > 0)
    {
        INFO("\t%d triangles", nr_ele_types[etype]);
    }
    if (nr_ele_types[++etype] > 0)
    {
        INFO("\t%d quads", nr_ele_types[etype]);
    }
    if (nr_ele_types[++etype] > 0)
    {
        INFO("\t%d tetrahedra", nr_ele_types[etype]);
    }
    if (nr_ele_types[++etype] > 0)
    {
        INFO("\t%d hexahedra", nr_ele_types[etype]);
    }
    if (nr_ele_types[++etype] > 0)
    {
        INFO("\t%d pyramids", nr_ele_types[etype]);
    }
    if (nr_ele_types[++etype] > 0)
    {
        INFO("\t%d prisms", nr_ele_types[etype]);
    }

    std::vector<std::string> const& vec_names (mesh->getProperties().getPropertyVectorNames());
    INFO("There are %d properties in the mesh:", vec_names.size());
    for (const auto& vec_name : vec_names)
    {
        if (auto const vec_bounds =
                MeshLib::MeshInformation::getValueBounds<int>(*mesh, vec_name))
        {
            INFO("\t%s: [%d, %d]", vec_name.c_str(), vec_bounds->first,
                 vec_bounds->second);
        }
        else if (auto const vec_bounds =
                     MeshLib::MeshInformation::getValueBounds<double>(*mesh,
                                                                      vec_name))
        {
            INFO("\t%s: [%g, %g]", vec_name.c_str(), vec_bounds->first,
                 vec_bounds->second);
        }
        else
        {
            INFO("\t%s: Could not get value bounds for property vector.",
                 vec_name.c_str());
        }
    }

    if (valid_arg.isSet()) {
        // MeshValidation outputs error messages
        // Remark: MeshValidation can modify the original mesh
        MeshLib::MeshValidation validation(*mesh);

        unsigned const n_holes (MeshLib::MeshValidation::detectHoles(*mesh));
        if (n_holes > 0)
        {
            INFO("%d hole(s) detected within the mesh", n_holes);
        }
        else
        {
            INFO ("No holes found within the mesh.");
        }
    }


    // *** write mesh in new format
    MeshLib::IO::writeMeshToFile(*mesh, ogs_mesh_arg.getValue());

    delete mesh;
}

