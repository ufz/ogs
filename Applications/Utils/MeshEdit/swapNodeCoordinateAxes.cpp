/**
 * @copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 */

#include <array>
#include <memory>
#include <string>

#include <tclap/CmdLine.h>

#include "Applications/ApplicationsLib/LogogSetup.h"

#include "BaseLib/BuildInfo.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"

#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/IO/writeMeshToFile.h"

static void swapNodeCoordinateAxes(MeshLib::Mesh &mesh, std::array<int, 3> const& new_axes_indices)
{
    double new_coords[3] = {};
    for (MeshLib::Node* node : mesh.getNodes())
    {
        for (int i=0; i<3; i++)
            new_coords[i] = (*node)[new_axes_indices[i]];
        for (int i=0; i<3; i++)
            (*node)[i] = new_coords[i];
    }
}

static bool parseNewOrder(std::string const& str_order, std::array<int, 3> &new_axes_indices)
{
    if (str_order.length()!=3)
    {
        ERR("Invalid argument for the new order. The argument should contain three characters.");
        return false;
    }

    for (std::size_t i=0; i<new_axes_indices.size(); i++)
        new_axes_indices[i] = -1;

    for (int i=0; i<3; i++)
    {
        if (str_order[i] == 'x')
            new_axes_indices[i] = 0;
        else if (str_order[i] == 'y')
            new_axes_indices[i] = 1;
        else if (str_order[i] == 'z')
            new_axes_indices[i] = 2;
        else {
            ERR("Invalid argument for the new order. The  given argument contains a character other than \"x\", \"y\", \"z\".");
            return false;
        }
    }

    bool isAxisSet[3] = {false};
    for (std::size_t i=0; i<new_axes_indices.size(); i++)
    {
        if (isAxisSet[new_axes_indices[i]])
        {
            ERR("Invalid argument for the new order. The argument contains some character used more than once.");
            return false;
        }
        isAxisSet[new_axes_indices[i]] = true;
    }

    return true;
}

int main(int argc, char *argv[])
{
    ApplicationsLib::LogogSetup logog_setup;

    TCLAP::CmdLine cmd("Swap node coordinate values", ' ', BaseLib::BuildInfo::git_describe);
    TCLAP::ValueArg<std::string> input_arg("i", "input-mesh-file","input mesh file",true,"","string");
    cmd.add( input_arg );
    TCLAP::ValueArg<std::string> output_arg("o", "output-mesh-file","output mesh file",true,"","string");
    cmd.add( output_arg );
    TCLAP::ValueArg<std::string> new_order_arg("n", "new-order", "the new order of swapped coordinate values "
                                                                 "(e.g. \"xzy\" for converting XYZ values to XZY values)",
                                               true, "", "string");
    cmd.add( new_order_arg );
    cmd.parse( argc, argv );

    const std::string str_order = new_order_arg.getValue();
    std::array<int, 3> new_order = {};
    if (!parseNewOrder(str_order, new_order))
        return EXIT_FAILURE;

    std::unique_ptr<MeshLib::Mesh> mesh(
        MeshLib::IO::readMeshFromFile(input_arg.getValue()));
    if (!mesh)
        return EXIT_FAILURE;

    if (mesh->getDimension() == 3)
    {
        WARN("Swapping coordinate values of 3D elements can result in incorrect node-ordering.");
    }

    INFO("Exchange node coordinates from xyz to %s", new_order_arg.getValue().data());
    swapNodeCoordinateAxes(*mesh, new_order);

    INFO("Save the new mesh into a file");
    MeshLib::IO::writeMeshToFile(*mesh, output_arg.getValue());

    return EXIT_SUCCESS;
}
