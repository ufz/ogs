/**
 * \file
 * \author Karsten Rink
 * \date   2011-12-19
 * \brief  Implementation of the generateMatPropsFromMatID tool.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <memory>

#include <tclap/CmdLine.h>

#include "Applications/ApplicationsLib/LogogSetup.h"

#include "BaseLib/BuildInfo.h"
#include "BaseLib/FileTools.h"

#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/IO/writeMeshToFile.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/Elements/Element.h"

int main (int argc, char* argv[])
{
    ApplicationsLib::LogogSetup logog_setup;

    TCLAP::CmdLine cmd(
        "Creates a new file for material properties and sets the material ids "
        "in the msh-file to 0\n\n"
        "OpenGeoSys-6 software, version " +
            BaseLib::BuildInfo::git_describe +
            ".\n"
            "Copyright (c) 2012-2018, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', BaseLib::BuildInfo::git_describe);

    TCLAP::ValueArg<std::string> mesh_arg("m",
                                              "mesh",
                                              "the mesh to open from a file",
                                              false,
                                              "",
                                              "filename for mesh input");
    cmd.add( mesh_arg );

    cmd.parse( argc, argv );

    // read mesh
    std::unique_ptr<MeshLib::Mesh> mesh(MeshLib::IO::readMeshFromFile(mesh_arg.getValue()));
    if (!mesh) {
        INFO("Could not read mesh from file \"%s\".", mesh_arg.getValue().c_str());
        return EXIT_FAILURE;
    }
    if (!mesh->getProperties().existsPropertyVector<int>("MaterialIDs"))
    {
        ERR("Mesh contains no int-property vector named \"MaterialIds\".");
        return EXIT_FAILURE;
    }
    auto materialIds = mesh->getProperties().getPropertyVector<int>("MaterialIDs");

    std::size_t const n_properties(materialIds->size());
    if (n_properties != mesh->getNumberOfElements()) {
        ERR("Size mismatch: number of elements (%u) != number of material "
            "properties (%u).", mesh->getNumberOfElements(), n_properties);
        return EXIT_FAILURE;
    }
    std::string const name = BaseLib::extractBaseNameWithoutExtension(mesh_arg.getValue());
    // create file
    std::string const new_matname(name + "_prop");
    std::ofstream out_prop( new_matname.c_str(), std::ios::out );
    if (out_prop.is_open())
    {
        for (std::size_t i=0; i<n_properties; ++i)
            out_prop << i << "\t" << (*materialIds)[i] << "\n";
        out_prop.close();
    }
    else
    {
        ERR("Could not create property \"%s\" file.", new_matname.c_str());
        return EXIT_FAILURE;
    }

    mesh->getProperties().removePropertyVector("MaterialIDs");

    std::string const new_mshname(name + "_new.vtu");
    INFO("Writing mesh to file \"%s\".", new_mshname.c_str());
    MeshLib::IO::writeMeshToFile(*mesh, new_mshname);

    INFO("New files \"%s\" and \"%s\" written.", new_mshname.c_str(), new_matname.c_str());

    return EXIT_SUCCESS;
}
