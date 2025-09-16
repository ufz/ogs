/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <tclap/CmdLine.h>
#include <vtkCleanUnstructuredGrid.h>
#include <vtkRemoveGhosts.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLPUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include "BaseLib/Logging.h"
#include "BaseLib/MPI.h"
#include "InfoLib/GitInfo.h"
#include "MeshLib/IO/writeMeshToFile.h"

int main(int argc, char* argv[])
{
    WARN(
        "Due to lack of functionality to handle VTU field data, this tool is "
        "replaced with a new tool, pvtu2vtu. Please use pvtu2vtu instead. If "
        "you use this tool, please make sure that the field data of the VTU "
        "file are not used in an OGS simulation.");

    TCLAP::CmdLine cmd(
        "Reads a VTK partitioned unstructured grid (*.pvtu), cleans the ghost "
        "information and saves the data as as a regular, connected mesh file."
        "\n\nOpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2025, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);

    TCLAP::ValueArg<std::string> output_arg(
        "o", "output", "Output (.vtu). The output mesh file", true, "",
        "OUTPUT_FILE");
    cmd.add(output_arg);

    TCLAP::ValueArg<std::string> input_arg(
        "i", "input", "Input (.pvtu). The partitioned input mesh file", true,
        "", "INPUT_FILE");
    cmd.add(input_arg);
    cmd.parse(argc, argv);

    BaseLib::MPI::Setup mpi_setup(argc, argv);

    vtkSmartPointer<vtkXMLPUnstructuredGridReader> reader =
        vtkSmartPointer<vtkXMLPUnstructuredGridReader>::New();
    reader->SetFileName(input_arg.getValue().c_str());

    vtkSmartPointer<vtkRemoveGhosts> ghosts =
        vtkSmartPointer<vtkRemoveGhosts>::New();
    ghosts->SetInputConnection(reader->GetOutputPort());

    vtkSmartPointer<vtkCleanUnstructuredGrid> clean =
        vtkSmartPointer<vtkCleanUnstructuredGrid>::New();
    clean->SetInputConnection(ghosts->GetOutputPort());

    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
        vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetInputConnection(clean->GetOutputPort());
    writer->SetFileName(output_arg.getValue().c_str());
    writer->Write();

    return EXIT_SUCCESS;
}
