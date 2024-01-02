/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <tclap/CmdLine.h>

#ifdef USE_PETSC
#include <mpi.h>
#endif

#include "Applications/FileIO/GocadIO/GocadAsciiReader.h"
#include "InfoLib/GitInfo.h"
#include "MeshLib/IO/VtkIO/VtuInterface.h"
#include "MeshLib/Mesh.h"

std::string getDelim(std::string const& str)
{
    std::size_t const bslash = str.find_first_of('\\');
    char const delim = (bslash == std::string::npos) ? '/' : '\\';
    return (str.back() == delim) ? "" : std::string(1, delim);
}

int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Reads Gocad ascii files (*.ts, *.pl, *.mx) and writes TSurf- and PLine"
        "data into one or more VTU unstructured grids.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2024, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);

    TCLAP::SwitchArg write_binary_arg(
        "b", "write-binary",
        "if set, OGS-Meshes will be written in binary format");
    cmd.add(write_binary_arg);

    TCLAP::SwitchArg export_surfaces_arg(
        "s", "surfaces-only",
        "if set, only TSurf datasets will be parsed from the input file");
    cmd.add(export_surfaces_arg);

    TCLAP::SwitchArg export_lines_arg(
        "l", "lines-only",
        "if set, only PLine datasets will be parsed from the input file");
    cmd.add(export_lines_arg);

    TCLAP::ValueArg<std::string> output_arg(
        "o", "output-dir", "output directory", true, "", "output dir");
    cmd.add(output_arg);

    TCLAP::ValueArg<std::string> input_arg(
        "i", "input-file",
        "Gocad triangular surfaces file (*.ts). Provide a file with unix file "
        "endings under unix. Use dos2unix to convert. ",
        true, "", "filename.ts");
    cmd.add(input_arg);

    cmd.parse(argc, argv);

#ifdef USE_PETSC
    MPI_Init(&argc, &argv);
#endif

    if (export_lines_arg.isSet() && export_surfaces_arg.isSet())
    {
        ERR("Both the 'lines-only'-flag and 'surfaces-only'-flag are set. Only "
            "one is allowed at a time.");
#ifdef USE_PETSC
        MPI_Finalize();
#endif
        return 2;
    }

    std::string const file_name(input_arg.getValue());

    FileIO::Gocad::DataType t(FileIO::Gocad::DataType::ALL);
    if (export_lines_arg.isSet())
    {
        t = FileIO::Gocad::DataType::PLINE;
    }
    if (export_surfaces_arg.isSet())
    {
        t = FileIO::Gocad::DataType::TSURF;
    }
    std::vector<std::unique_ptr<MeshLib::Mesh>> meshes;
    if (!FileIO::Gocad::GocadAsciiReader::readFile(file_name, meshes, t))
    {
        ERR("Error reading file.");
#ifdef USE_PETSC
        MPI_Finalize();
#endif
        return 1;
    }
    INFO("{:d} meshes found.", meshes.size());
    std::string const dir = output_arg.getValue();
    bool const write_binary = write_binary_arg.getValue();
    std::string const delim = getDelim(dir);
    for (auto& mesh : meshes)
    {
        if (mesh == nullptr)
        {
            continue;
        }
        INFO("Writing mesh \"{:s}\"", mesh->getName());
        int data_mode = (write_binary) ? 2 : 0;
        bool compressed = (write_binary);
        MeshLib::IO::VtuInterface vtu(mesh.get(), data_mode, compressed);
        vtu.writeToFile(dir + delim + mesh->getName() + ".vtu");
    }
#ifdef USE_PETSC
    MPI_Finalize();
#endif
    return 0;
}
