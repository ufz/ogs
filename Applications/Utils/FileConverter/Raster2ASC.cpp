/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <tclap/CmdLine.h>

#include "BaseLib/Logging.h"
#include "BaseLib/MPI.h"
#include "BaseLib/TCLAPArguments.h"
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
            "Copyright (c) 2012-2025, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);

    TCLAP::ValueArg<std::string> output_arg(
        "o", "output", "Output (.asc). Filename for output raster", true, "",
        "OUTPUT_FILE");
    cmd.add(output_arg);

    TCLAP::ValueArg<std::string> input_arg(
        "i", "input", "Input (.grd | .xyz). Filename for input raster", true,
        "", "INPUT_FILE");
    cmd.add(input_arg);

    auto log_level_arg = BaseLib::makeLogLevelArg();
    cmd.add(log_level_arg);
    cmd.parse(argc, argv);

    BaseLib::MPI::Setup mpi_setup(argc, argv);
    BaseLib::initOGSLogger(log_level_arg.getValue());

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
    return EXIT_SUCCESS;
}
