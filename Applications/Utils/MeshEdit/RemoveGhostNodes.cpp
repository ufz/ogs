/**
 * \file
 *
 * @copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <tclap/CmdLine.h>

#include "Applications/ApplicationsLib/LogogSetup.h"
#include "InfoLib/GitInfo.h"

#include "MeshLib/IO/writeMeshToFile.h"
#include "vtkCleanUnstructuredGrid.h"

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLPUnstructuredGridReader.h>
#include <vtkRemoveGhosts.h>

int main (int argc, char* argv[])
{
    ApplicationsLib::LogogSetup logog_setup;
/*
    TCLAP::CmdLine cmd(
        "Reads a 3D unstructured mesh and samples it onto a structured grid of "
        "the same extent. Cell properties are mapped onto the grid (sampled at "
        "the centre-points of each cube), node properties are ignored. Note, "
        "that a large cube size may result in an undersampling of the original "
        "mesh structure.\nCube sizes are defines by x/y/z-parameters. For "
        "equilateral cubes, only the x-parameter needs to be set.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2020, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);

    TCLAP::ValueArg<double> z_arg("z", "cellsize-z",
                                  "edge length of cubes in z-direction (depth)",
                                  false, 1000, "floating point number");
    cmd.add(z_arg);

    TCLAP::ValueArg<double> y_arg(
        "y", "cellsize-y", "edge length of cubes in y-direction (latitude)",
        false, 1000, "floating point number");
    cmd.add(y_arg);

    TCLAP::ValueArg<double> x_arg(
        "x", "cellsize-x",
        "edge length of cubes in x-direction (longitude) or all directions, if "
        "y and z are not set",
        true, 1000, "floating point number");
    cmd.add(x_arg);

    TCLAP::ValueArg<std::string> output_arg(
        "o", "output", "the output grid (*.vtu)", true, "", "output.vtu");
    cmd.add(output_arg);

    TCLAP::ValueArg<std::string> input_arg("i", "input",
                                           "the 3D input mesh (*.vtu, *.msh)",
                                           true, "", "input.vtu");
    cmd.add(input_arg);
    cmd.parse(argc, argv);

    if ((y_arg.isSet() && !z_arg.isSet()) ||
        ((!y_arg.isSet() && z_arg.isSet())))
    {
        ERR("For equilateral cubes, only x needs to be set. For unequal "
            "cuboids, all three edge lengths (x/y/z) need to be specified.")
        return -1;
    }

    double const x_size = x_arg.getValue();
    double const y_size = (y_arg.isSet()) ? y_arg.getValue() : x_arg.getValue();
    double const z_size = (z_arg.isSet()) ? z_arg.getValue() : x_arg.getValue();
    std::array<double, 3> const cellsize = { x_size, y_size, z_size };
*/
    std::string input_file = "c:/Projects/RemoveGhostNodes/Mesh3D.pvtu";
    vtkSmartPointer<vtkXMLPUnstructuredGridReader> reader =
        vtkSmartPointer<vtkXMLPUnstructuredGridReader>::New();
    reader->SetFileName(input_file.c_str());
    //reader->SetFileName(input_arg.getValue().c_str());
    reader->Update();
    vtkSmartPointer<vtkUnstructuredGrid> mesh = reader->GetOutput();

    vtkSmartPointer<vtkRemoveGhosts> ghosts =
        vtkSmartPointer<vtkRemoveGhosts>::New();
    ghosts->SetInputConnection(reader->GetOutputPort());


    //if (MeshLib::IO::writeMeshToFile(*grid, output_arg.getValue()) != 0)
    //    return EXIT_FAILURE;

    return EXIT_SUCCESS;
}
