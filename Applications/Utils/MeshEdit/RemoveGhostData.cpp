/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
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

#include "InfoLib/GitInfo.h"
#include "MeshLib/IO/writeMeshToFile.h"

int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Reads a VTK partitioned unstructured grid (*.pvtu), cleans the ghost "
        "information and saves the data as as a regular, connected mesh file."
        "\n\nOpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2021, OpenGeoSys Community "
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
