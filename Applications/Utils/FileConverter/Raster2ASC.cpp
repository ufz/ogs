/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <tclap/CmdLine.h>

#ifdef USE_PETSC
#include <mpi.h>
#endif

#include "GeoLib/IO/AsciiRasterInterface.h"
#include "GeoLib/Raster.h"
#include "InfoLib/GitInfo.h"

int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Converts ascii raster files (e.g. Surfer *.grd files or *.xyz files) "
        "into ASC raster files.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2024, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);

    TCLAP::ValueArg<std::string> output_arg(
        "o", "output", "filename for output raster", true, "", "output file");
    cmd.add(output_arg);

    TCLAP::ValueArg<std::string> input_arg(
        "i", "input", "filename for input raster", true, "", "input file");
    cmd.add(input_arg);

    cmd.parse(argc, argv);

#ifdef USE_PETSC
    MPI_Init(&argc, &argv);
#endif

    std::unique_ptr<GeoLib::Raster> raster(
        FileIO::AsciiRasterInterface::readRaster(input_arg.getValue()));

    if (raster == nullptr)
    {
        ERR("Couldn't read input raster file.");
        return EXIT_FAILURE;
    }

    std::string output_name = output_arg.getValue();
    if (output_name.substr(output_name.length() - 4, 4) != ".asc")
    {
        output_name = output_name.append(".asc");
    }

    FileIO::AsciiRasterInterface::writeRasterAsASC(*raster, output_name);
#ifdef USE_PETSC
    MPI_Finalize();
#endif
    return EXIT_SUCCESS;
}
