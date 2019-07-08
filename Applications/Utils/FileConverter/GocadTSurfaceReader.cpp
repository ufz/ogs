/**
 *
 * @copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <tclap/CmdLine.h>

#include "BaseLib/BuildInfo.h"
#include "Applications/ApplicationsLib/LogogSetup.h"
#include "Applications/FileIO/GocadIO/GocadTSurfaceReader.h"

int main(int argc, char* argv[])
{
    ApplicationsLib::LogogSetup logog_setup;

    TCLAP::CmdLine cmd(
        "Reads a Gocad triangular surfaces file (*.ts) and writes the "
        "data into one or more VTU unstructured grids.\n\n"
        "OpenGeoSys-6 software, version " +
            BaseLib::BuildInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2019, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', BaseLib::BuildInfo::ogs_version);

    TCLAP::ValueArg<std::string> input_arg(
        "i", "input-file", "Gocad triangular surfaces file (*.ts)", true, "",
        "filename.ts");
    cmd.add(input_arg);

    TCLAP::ValueArg<std::string> output_arg(
        "o", "output-dir", "output directory", true, "",
        "output dir");
    cmd.add(output_arg);

    TCLAP::SwitchArg write_binary_arg(
        "b", "write-binary",
        "if set, OGS-Meshes will be written in binary format");
    cmd.add(write_binary_arg);

    cmd.parse(argc, argv);

    FileIO::Gocad::GocadTSurfaceReader gcts(input_arg.getValue());
    gcts.readFile();
    gcts.writeData(output_arg.getValue(), write_binary_arg.getValue());

    return 0;
}
