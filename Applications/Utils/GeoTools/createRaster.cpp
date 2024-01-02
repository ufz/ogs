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

#include <memory>

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
            "Copyright (c) 2012-2024, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);

    TCLAP::ValueArg<std::string> output_arg("o", "output",
                                            "Name of the output raster (*.asc)",
                                            true, "", "output file name");
    cmd.add(output_arg);
    TCLAP::ValueArg<std::size_t> n_rows("r", "n_rows", "number of rows", false,
                                        1000, "positive integer value");
    cmd.add(n_rows);
    TCLAP::ValueArg<std::size_t> n_cols("c",
                                        "n_cols",
                                        "number of columns",
                                        false,
                                        1000,
                                        "positive integer value");
    cmd.add(n_cols);
    TCLAP::ValueArg<double> cell_size("s", "cell_size", "cell size", false,
                                      10.0, "double value");
    cmd.add(cell_size);
    TCLAP::ValueArg<double> ll_y_arg(
        "",
        "ll_y",
        "y coordinate of lower left point of axis aligned rectangular region",
        false,
        0,
        "double value");
    cmd.add(ll_y_arg);
    TCLAP::ValueArg<double> ll_x_arg(
        "",
        "ll_x",
        "x coordinate of lower left point of axis aligned rectangular region",
        false,
        0,
        "double value");
    cmd.add(ll_x_arg);

    cmd.parse(argc, argv);

#ifdef USE_PETSC
    MPI_Init(&argc, &argv);
#endif

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

#ifdef USE_PETSC
    MPI_Finalize();
#endif
    return EXIT_SUCCESS;
}
