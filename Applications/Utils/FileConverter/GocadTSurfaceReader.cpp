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
#include "MeshLib/Mesh.h"
#include "MeshLib/IO/VtkIO/VtuInterface.h"
#include "Applications/ApplicationsLib/LogogSetup.h"
#include "Applications/FileIO/GocadIO/GocadTSurfaceReader.h"

std::string getDelim(std::string const& str)
{
    std::size_t const bslash = str.find_first_of('\\');
    char const delim = (bslash == std::string::npos) ? '/' : '\\';
    return (str.back() == delim) ? "" : std::string(1, delim);
}

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

    std::string const file_name (input_arg.getValue());
    FileIO::Gocad::GocadTSurfaceReader gcts;
    std::vector<std::unique_ptr<MeshLib::Mesh>> meshes;
    if (!gcts.readFile(file_name, meshes))
    {
        ERR("Error reading file.");
        return 1;
    }
    std::string const dir = output_arg.getValue();
    bool const write_binary = write_binary_arg.getValue();
    std::string const delim = getDelim(dir);
    for (auto& mesh : meshes)
    {
        INFO("Writing mesh \"%s\"", mesh->getName().c_str());
        int data_mode = (write_binary) ? 2 : 0;
        bool compressed = (write_binary);
        MeshLib::IO::VtuInterface vtu(mesh.get(), data_mode, compressed);
        vtu.writeToFile(dir + delim + mesh->getName() + ".vtu");
    }
    return 0;
}
