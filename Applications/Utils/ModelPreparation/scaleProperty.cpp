/**
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <algorithm>
#include <memory>
#include <cmath>
#include <numeric>

#include <tclap/CmdLine.h>

#include "Applications/ApplicationsLib/LogogSetup.h"

#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/IO/writeMeshToFile.h"

#include "MeshLib/Mesh.h"

int main(int argc, char* argv[])
{
    ApplicationsLib::LogogSetup logo_setup;

    TCLAP::CmdLine cmd(
        "Scales a property of a mesh",
        ' ',
        "0.1");

    TCLAP::ValueArg<std::string> out_mesh_arg(
        "o",
        "out-mesh",
        "the mesh is stored to a file of this name",
        true,
        "",
        "filename for mesh output");
    cmd.add(out_mesh_arg);

    TCLAP::ValueArg<std::string> property_arg(
        "p",
        "property-name",
        "the name of the property the values are stored for",
        true,
        "",
        "property name as string");
    cmd.add(property_arg);

    TCLAP::ValueArg<std::string> mesh_arg("m",
                                          "mesh",
                                          "the mesh is read from this file",
                                          true,
                                          "",
                                          "file name");
    cmd.add(mesh_arg);

    std::vector<std::string> allowed_units{ "mm/a", "mm/month", "m/s" };
    TCLAP::ValuesConstraint<std::string> allowed_units_constraints{
        allowed_units};
    TCLAP::ValueArg<std::string> unit_arg("u",
                                          "input-unit",
                                          "input unit of the data",
                                          true,
                                          "m/s",
                                          &allowed_units_constraints);
    cmd.add(unit_arg);

    cmd.parse(argc, argv);

    std::unique_ptr<MeshLib::Mesh> mesh(
        MeshLib::IO::readMeshFromFile(mesh_arg.getValue()));

    double scale(1.0);
    if (unit_arg.getValue() == "m/s")
    {
        scale = 1.0;
    }
    else if (unit_arg.getValue() == "mm/a")
    {
        scale = 1e-3 / (365.25 * 86400);
    }
    else if (unit_arg.getValue() == "mm/month")
    {
        scale = 1e-3 * (12.0 / (365.25 * 86400));
    }

    MeshLib::scaleMeshPropertyVector(*mesh, property_arg.getValue(), scale);

    MeshLib::IO::writeMeshToFile(*mesh, out_mesh_arg.getValue());

    return EXIT_SUCCESS;
}
