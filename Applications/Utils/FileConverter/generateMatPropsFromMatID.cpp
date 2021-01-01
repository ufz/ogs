/**
 * \file
 * \author Karsten Rink
 * \date   2011-12-19
 * \brief  Implementation of the generateMatPropsFromMatID tool.
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <memory>

#include <tclap/CmdLine.h>

#include "InfoLib/GitInfo.h"
#include "BaseLib/FileTools.h"

#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/IO/writeMeshToFile.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/Elements/Element.h"

int main (int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Creates a new file for material properties and sets the material ids "
        "in the msh-file to 0\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2021, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);

    TCLAP::ValueArg<std::string> mesh_arg("m",
                                              "mesh",
                                              "the mesh to open from a file",
                                              false,
                                              "",
                                              "filename for mesh input");
    cmd.add( mesh_arg );

    cmd.parse( argc, argv );

    // read mesh
    std::unique_ptr<MeshLib::Mesh> mesh(
        MeshLib::IO::readMeshFromFile(mesh_arg.getValue()));

    if (!mesh)
    {
        INFO("Could not read mesh from file '{:s}'.", mesh_arg.getValue());
        return EXIT_FAILURE;
    }

    auto const materialIds = materialIDs(*mesh);
    if (!materialIds)
    {
        OGS_FATAL("Mesh contains no int-property vector named 'MaterialIDs'.");
    }

    std::size_t const n_properties(materialIds->size());
    assert(n_properties != mesh->getNumberOfElements());

    std::string const name =
        BaseLib::extractBaseNameWithoutExtension(mesh_arg.getValue());
    // create file
    std::string const new_matname(name + "_prop");
    std::ofstream out_prop(new_matname.c_str(), std::ios::out);
    if (out_prop.is_open())
    {
        for (std::size_t i = 0; i < n_properties; ++i)
        {
            out_prop << i << "\t" << (*materialIds)[i] << "\n";
        }
        out_prop.close();
    }
    else
    {
        ERR("Could not create property '{:s}' file.", new_matname);
        return EXIT_FAILURE;
    }

    mesh->getProperties().removePropertyVector("MaterialIDs");

    std::string const new_mshname(name + "_new.vtu");
    INFO("Writing mesh to file '{:s}'.", new_mshname);
    MeshLib::IO::writeMeshToFile(*mesh, new_mshname);

    INFO("New files '{:s}' and '{:s}' written.", new_mshname, new_matname);

    return EXIT_SUCCESS;
}
