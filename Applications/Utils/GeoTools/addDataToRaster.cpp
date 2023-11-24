/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <tclap/CmdLine.h>

#ifdef USE_PETSC
#include <mpi.h>
#endif

#include <algorithm>
#include <boost/math/constants/constants.hpp>
#include <cmath>
#include <memory>
#include <numeric>

#include "GeoLib/AABB.h"
#include "GeoLib/IO/AsciiRasterInterface.h"
#include "GeoLib/Point.h"
#include "GeoLib/Raster.h"
#include "InfoLib/GitInfo.h"

double compute2DGaussBellCurveValues(GeoLib::Point const& point,
                                     GeoLib::AABB const& aabb)
{
    auto const sigma_x = (aabb.getMaxPoint() - aabb.getMinPoint())[0] / 3;
    auto const sigma_y = (aabb.getMaxPoint() - aabb.getMinPoint())[1] / 3;

    auto const mid_point = (aabb.getMaxPoint() + aabb.getMinPoint()) / 2;

    return std::exp(
        -0.5 * std::pow((point[0] - mid_point[0]), 2) / std::pow(sigma_x, 2) -
        0.5 * std::pow((point[1] - mid_point[1]), 2) / std::pow(sigma_y, 2));
}

double computeSinXSinY(GeoLib::Point const& point, GeoLib::AABB const& aabb)
{
    auto const aabb_size = aabb.getMaxPoint() - aabb.getMinPoint();
    auto const offset = aabb.getMinPoint();

    return std::sin((point[0] - offset[0]) / aabb_size[0] *
                    boost::math::double_constants::pi) *
           std::sin((point[1] - offset[1]) / aabb_size[1] *
                    boost::math::double_constants::pi);
}

double computeStepFunction(GeoLib::Point const& point, GeoLib::AABB const& aabb)
{
    auto const aabb_size = aabb.getMaxPoint() - aabb.getMinPoint();
    auto const min = aabb.getMinPoint();

    if ((point[0] - min[0]) / aabb_size[0] < 0.5 &&
        (point[1] - min[1]) / aabb_size[1] < 0.5)
    {
        return 1;
    }
    if ((point[0] - min[0]) / aabb_size[0] >= 0.5 &&
        (point[1] - min[1]) / aabb_size[1] >= 0.5)
    {
        return 1;
    }
    return 0;
}

int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Add values to raster.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2023, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);

    TCLAP::ValueArg<std::string> out_raster_arg(
        "o",
        "output_raster",
        "the output raster is stored to a file of this name",
        true,
        "",
        "filename for raster output");
    cmd.add(out_raster_arg);

    TCLAP::ValueArg<double> scaling_arg(
        "",
        "scaling_value",
        "value the function sin(x pi) sin(y pi) will be scaled with",
        false,
        1,
        "double value");
    cmd.add(scaling_arg);

    TCLAP::ValueArg<double> offset_arg(
        "",
        "offset_value",
        "constant added to the function 'scaling * sin(x pi) * sin(y pi)'",
        false,
        0,
        "double value");
    cmd.add(offset_arg);

    TCLAP::ValueArg<double> ll_x_arg(
        "",
        "ll_x",
        "x coordinate of lower left point of axis aligned rectangular region",
        false,
        0,
        "double value");
    cmd.add(ll_x_arg);
    TCLAP::ValueArg<double> ll_y_arg(
        "",
        "ll_y",
        "y coordinate of lower left point of axis aligned rectangular region",
        false,
        0,
        "double value");
    cmd.add(ll_y_arg);
    TCLAP::ValueArg<double> ur_x_arg("",
                                     "ur_x",
                                     "x coordinate of the upper right point of "
                                     "axis aligned rectangular region",
                                     false,
                                     0,
                                     "double value");
    cmd.add(ur_x_arg);
    TCLAP::ValueArg<double> ur_y_arg("",
                                     "ur_y",
                                     "y coordinate of the upper right point of "
                                     "axis aligned rectangular region",
                                     false,
                                     0,
                                     "double value");

    cmd.add(ur_y_arg);
    std::vector<std::string> allowed_functions_vector{"sinxsiny", "exp",
                                                      "step"};
    TCLAP::ValuesConstraint<std::string> allowed_functions(
        allowed_functions_vector);
    TCLAP::ValueArg<std::string> function_arg(
        "f", "function", "Name of the function used to modify the raster", true,
        "", &allowed_functions);
    cmd.add(function_arg);
    TCLAP::ValueArg<std::string> input_arg("i", "input",
                                           "Name of the input raster (*.asc)",
                                           true, "", "input file name");
    cmd.add(input_arg);

    cmd.parse(argc, argv);

#ifdef USE_PETSC
    MPI_Init(&argc, &argv);
#endif

    std::array input_points = {
        GeoLib::Point{{ll_x_arg.getValue(), ll_y_arg.getValue(), 0}},
        GeoLib::Point{{ur_x_arg.getValue(), ur_y_arg.getValue(), 0}}};
    GeoLib::AABB const aabb{std::begin(input_points), std::end(input_points)};

    auto const s = scaling_arg.getValue();
    auto const offset = offset_arg.getValue();

    std::unique_ptr<GeoLib::Raster> const raster(
        FileIO::AsciiRasterInterface::getRasterFromASCFile(
            input_arg.getValue()));
    auto const& header = raster->getHeader();
    auto const& origin = header.origin;

    auto function_selector = [](std::string const& function_string)
        -> std::function<double(GeoLib::Point const& p,
                                GeoLib::AABB const& aabb)>
    {
        if (function_string == "sinxsiny")
        {
            return computeSinXSinY;
        }
        if (function_string == "exp")
        {
            return compute2DGaussBellCurveValues;
        }
        if (function_string == "step")
        {
            return computeStepFunction;
        }
        OGS_FATAL("Function '{}' isn't implemented.", function_string);
    };
    std::function<double(GeoLib::Point const& p, GeoLib::AABB const& aabb)>
        computeFunctionValue = function_selector(function_arg.getValue());

    for (std::size_t r = 0; r < header.n_rows; r++)
    {
        for (std::size_t c = 0; c < header.n_cols; c++)
        {
            GeoLib::Point const p{{origin[0] + header.cell_size * (c + 0.5),
                                   origin[1] + header.cell_size * (r + 0.5),
                                   0.0}};
            if (!aabb.containsPoint(p, std::numeric_limits<double>::epsilon()))
            {
                continue;
            }

            (*raster)(r, c) += offset + s * computeFunctionValue(p, aabb);
        }
    }

    FileIO::AsciiRasterInterface::writeRasterAsASC(*raster,
                                                   out_raster_arg.getValue());
#ifdef USE_PETSC
    MPI_Finalize();
#endif
    return EXIT_SUCCESS;
}
