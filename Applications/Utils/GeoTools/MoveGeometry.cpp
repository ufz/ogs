// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <tclap/CmdLine.h>

#include "BaseLib/Logging.h"
#include "BaseLib/MPI.h"
#include "BaseLib/TCLAPArguments.h"
#include "GeoLib/GEOObjects.h"
#include "GeoLib/IO/XmlIO/Boost/BoostXmlGmlInterface.h"
#include "InfoLib/GitInfo.h"

int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Moves the points of a geometry by a given displacement vector\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2026, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);
    TCLAP::ValueArg<double> z_arg("z", "z", "displacement in z direction",
                                  false, 0.0, "Z-DISPLACEMENT");
    cmd.add(z_arg);
    TCLAP::ValueArg<double> y_arg("y", "y", "displacement in y direction",
                                  false, 0.0, "Y-DISPLACEMENT");
    cmd.add(y_arg);
    TCLAP::ValueArg<double> x_arg("x", "x", "displacement in x direction",
                                  false, 0.0, "X-DISPLACEMENT");
    cmd.add(x_arg);
    TCLAP::ValueArg<std::string> geo_output_arg(
        "o", "output", "Output (.gml) geometry file", true, "", "OUTPUT_FILE");
    cmd.add(geo_output_arg);
    TCLAP::ValueArg<std::string> geo_input_arg(
        "i", "input", "Input (.gml) geometry file", true, "", "INPUT_FILE");
    cmd.add(geo_input_arg);
    auto log_level_arg = BaseLib::makeLogLevelArg();
    cmd.add(log_level_arg);
    cmd.parse(argc, argv);

    BaseLib::MPI::Setup mpi_setup(argc, argv);
    BaseLib::initOGSLogger(log_level_arg.getValue());

    GeoLib::GEOObjects geo_objects;
    GeoLib::IO::BoostXmlGmlInterface xml(geo_objects);
    try
    {
        if (!xml.readFile(geo_input_arg.getValue()))
        {
            return EXIT_FAILURE;
        }
    }
    catch (std::runtime_error const& err)
    {
        ERR("Failed to read file `{:s}'.", geo_input_arg.getValue());
        ERR("{:s}", err.what());
        return EXIT_FAILURE;
    }

    Eigen::Vector3d displacement = Eigen::Vector3d::Zero();
    if (x_arg.isSet())
    {
        displacement[0] = x_arg.getValue();
    }
    if (y_arg.isSet())
    {
        displacement[1] = y_arg.getValue();
    }
    if (z_arg.isSet())
    {
        displacement[2] = z_arg.getValue();
    }

    auto const geo_name = geo_objects.getGeometryNames()[0];

    std::vector<GeoLib::Point*> const* point_vec =
        geo_objects.getPointVec(geo_name);
    std::size_t const n_points = point_vec->size();
    for (std::size_t i = 0; i < n_points; ++i)
    {
        for (std::size_t c = 0; c < 3; ++c)
        {
            (*(*point_vec)[i])[c] += displacement[c];
        }
    }

    xml.export_name = geo_name;
    BaseLib::IO::writeStringToFile(xml.writeToString(),
                                   geo_output_arg.getValue());

    return EXIT_SUCCESS;
}
