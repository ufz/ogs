/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
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

#ifdef USE_PETSC
#include <mpi.h>
#endif

#include "BaseLib/Logging.h"
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
            "Copyright (c) 2012-2024, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);

    TCLAP::ValueArg<std::string> output_arg(
        "o", "output", "the output mesh (*.vtu)", true, "", "output.vtu");
    cmd.add(output_arg);

    TCLAP::ValueArg<std::string> input_arg(
        "i", "input", "the partitioned input mesh (*.pvtu)", true, "",
        "input.pvtu");
    cmd.add(input_arg);
    cmd.parse(argc, argv);

#ifdef USE_PETSC
    MPI_Init(&argc, &argv);
#endif

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

#ifdef USE_PETSC
    MPI_Finalize();
#endif
    return EXIT_SUCCESS;
}
