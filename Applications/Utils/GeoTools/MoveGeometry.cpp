/**
 * \file
 * \author Karsten Rink
 * \date   2015-08-04
 * \brief  A small tool to translate geometries.
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

// ThirdParty
#include <tclap/CmdLine.h>

#ifdef USE_PETSC
#include <mpi.h>
#endif

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
            "Copyright (c) 2012-2024, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);
    TCLAP::ValueArg<double> z_arg("z", "z", "displacement in z direction",
                                  false, 0.0, "z-displacement");
    cmd.add(z_arg);
    TCLAP::ValueArg<double> y_arg("y", "y", "displacement in y direction",
                                  false, 0.0, "y-displacement");
    cmd.add(y_arg);
    TCLAP::ValueArg<double> x_arg("x", "x", "displacement in x direction",
                                  false, 0.0, "x-displacement");
    cmd.add(x_arg);
    TCLAP::ValueArg<std::string> geo_output_arg(
        "o", "output", "output geometry file (*.gml)", true, "", "output file");
    cmd.add(geo_output_arg);
    TCLAP::ValueArg<std::string> geo_input_arg(
        "i", "input", "input geometry file (*.gml)", true, "", "input file");
    cmd.add(geo_input_arg);
    cmd.parse(argc, argv);

#ifdef USE_PETSC
    MPI_Init(&argc, &argv);
#endif

    GeoLib::GEOObjects geo_objects;
    GeoLib::IO::BoostXmlGmlInterface xml(geo_objects);
    try
    {
        if (!xml.readFile(geo_input_arg.getValue()))
        {
#ifdef USE_PETSC
            MPI_Finalize();
#endif
            return EXIT_FAILURE;
        }
    }
    catch (std::runtime_error const& err)
    {
        ERR("Failed to read file `{:s}'.", geo_input_arg.getValue());
        ERR("{:s}", err.what());
#ifdef USE_PETSC
        MPI_Finalize();
#endif
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

#ifdef USE_PETSC
    MPI_Finalize();
#endif
    return EXIT_SUCCESS;
}
