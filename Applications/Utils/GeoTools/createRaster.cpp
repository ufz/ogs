// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <tclap/CmdLine.h>

#include <memory>

#include "BaseLib/Logging.h"
#include "BaseLib/MPI.h"
#include "BaseLib/TCLAPArguments.h"
#include "BaseLib/TCLAPOutput.h"
#include "GeoLib/AABB.h"
#include "GeoLib/IO/AsciiRasterInterface.h"
#include "GeoLib/Point.h"
#include "GeoLib/Raster.h"
#include "InfoLib/GitInfo.h"

int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Create a raster of specified size at specified origin where every "
        "pixel has the value zero.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2025, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);
    BaseLib::TCLAPOutput tclapOutput;
    cmd.setOutput(&tclapOutput);

    TCLAP::ValueArg<std::string> output_arg("o", "output",
                                            "Output (.asc). Name of the output"
                                            "raster file",
                                            true, "", "OUTPUT_FILE");
    cmd.add(output_arg);
    TCLAP::ValueArg<std::size_t> n_rows("r", "n_rows", "number of rows", false,
                                        1000, "NUM_ROWS");
    cmd.add(n_rows);
    TCLAP::ValueArg<std::size_t> n_cols("c", "n_cols", "number of columns",
                                        false, 1000, "NUM_COLS");
    cmd.add(n_cols);
    TCLAP::ValueArg<double> cell_size("s", "cell_size",
                                      "cell size, "
                                      "(min = 0)",
                                      false, 10.0, "CELL_SIZE");
    cmd.add(cell_size);
    TCLAP::ValueArg<double> ll_y_arg(
        "",
        "ll_y",
        "y coordinate of lower left point of axis aligned rectangular region",
        false,
        0,
        "LL_Y");
    cmd.add(ll_y_arg);
    TCLAP::ValueArg<double> ll_x_arg(
        "",
        "ll_x",
        "x coordinate of lower left point of axis aligned rectangular region",
        false,
        0,
        "LL_X");
    cmd.add(ll_x_arg);

    auto log_level_arg = BaseLib::makeLogLevelArg();
    cmd.add(log_level_arg);
    cmd.parse(argc, argv);

    BaseLib::MPI::Setup mpi_setup(argc, argv);
    BaseLib::initOGSLogger(log_level_arg.getValue());

    GeoLib::RasterHeader header{
        n_cols.getValue(),
        n_rows.getValue(),
        0,
        GeoLib::Point{{ll_x_arg.getValue(), ll_y_arg.getValue(), 0}},
        cell_size.getValue(),
        -9999};
    std::vector<double> raster_data(header.n_cols * header.n_rows, 0.0);
    GeoLib::Raster const raster{header, raster_data.begin(), raster_data.end()};

    FileIO::AsciiRasterInterface::writeRasterAsASC(raster,
                                                   output_arg.getValue());

    return EXIT_SUCCESS;
}
