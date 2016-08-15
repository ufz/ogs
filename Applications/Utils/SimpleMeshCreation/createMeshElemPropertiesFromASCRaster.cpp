/**
 * \brief  Implementation of the createMeshElemPropertiesFromASCRaster tool.
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <algorithm>
#include <memory>
#include <cmath>
#include <numeric>

#include <tclap/CmdLine.h>

#include "Applications/ApplicationsLib/LogogSetup.h"

#include "BaseLib/quicksort.h"
#include "BaseLib/FileTools.h"

#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/IO/writeMeshToFile.h"
#include "Applications/FileIO/AsciiRasterInterface.h"

#include "GeoLib/Raster.h"

#include "MathLib/MathTools.h"

#include "MeshLib/MeshGenerators/RasterToMesh.h"
#include "MeshLib/MeshGenerators/VtkMeshConverter.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshEditing/Mesh2MeshPropertyInterpolation.h"
#include "MeshLib/MeshEnums.h"

void scaleMeshPropertyVector(MeshLib::Mesh & mesh,
                             std::string const& property_name,
                             double factor)
{
    boost::optional<MeshLib::PropertyVector<double> &> pv(
        mesh.getProperties().getPropertyVector<double>(property_name));

    for (auto & v : *pv)
        v *= factor;
}

int main(int argc, char* argv[])
{
    ApplicationsLib::LogogSetup logo_setup;

    TCLAP::CmdLine cmd(
        "Generates properties for mesh elements of an input mesh deploying a "
        "ASC raster file",
        ' ',
        "0.1");

    TCLAP::ValueArg<std::string> out_mesh_arg(
        "o",
        "out-mesh",
        "the mesh is stored to a file of this name",
        false,
        "",
        "filename for mesh output");
    cmd.add(out_mesh_arg);

    TCLAP::ValueArg<bool> refinement_raster_output_arg(
        "",
        "output-refined-raster",
        "write refined raster to a new ASC file",
        false,
        false,
        "0");
    cmd.add(refinement_raster_output_arg);

    TCLAP::ValueArg<unsigned> refinement_arg(
        "r",
        "refine",
        "refinement factor that raises the resolution of the raster data",
        false,
        1,
        "factor (default = 1)");
    cmd.add(refinement_arg);

    TCLAP::ValueArg<std::string> raster_arg("",
                                            "raster-file",
                                            "the name of the ASC raster file",
                                            true,
                                            "",
                                            "file name");
    cmd.add(raster_arg);

    TCLAP::ValueArg<std::string> property_arg(
        "p",
        "property-name",
        "the name of the property the values are stored for",
        true,
        "",
        "property name as string");
    cmd.add(property_arg);

    TCLAP::ValueArg<std::string> mesh_arg("m",
                                          "mesh",
                                          "the mesh is read from this file",
                                          true,
                                          "",
                                          "file name");
    cmd.add(mesh_arg);

    std::vector<std::string> allowed_units{ "mm/a", "mm/month", "m/s" };
    TCLAP::ValuesConstraint<std::string> allowed_units_constraints{
        allowed_units};
    TCLAP::ValueArg<std::string> unit_arg("u",
                                          "input-unit",
                                          "input unit of the data",
                                          true,
                                          "m/s",
                                          &allowed_units_constraints);
    cmd.add(unit_arg);

    cmd.parse(argc, argv);

    // read mesh
    std::unique_ptr<MeshLib::Mesh> dest_mesh(
        MeshLib::IO::readMeshFromFile(mesh_arg.getValue()));

    // read raster and if required manipulate it
    auto raster = std::unique_ptr<GeoLib::Raster>(
        FileIO::AsciiRasterInterface::getRasterFromASCFile(
            raster_arg.getValue()));
    GeoLib::RasterHeader header(raster->getHeader());
    if (refinement_arg.getValue() > 1)
    {
        raster->refineRaster(refinement_arg.getValue());
        if (refinement_raster_output_arg.getValue())
        {
            // write new asc file
            std::string new_raster_fname(
                BaseLib::dropFileExtension(raster_arg.getValue()));
            new_raster_fname += "-" + std::to_string(header.n_rows) + "x" +
                                std::to_string(header.n_cols) + ".asc";
            FileIO::AsciiRasterInterface::writeRasterAsASC(*raster,
                                                           new_raster_fname);
        }
    }

    // put raster data in a std::vector
    GeoLib::Raster::const_iterator raster_it(raster->begin());
    std::size_t size(header.n_cols * header.n_rows);
    std::vector<double> src_properties(size);
    for (unsigned row(0); row < header.n_rows; row++)
    {
        for (unsigned col(0); col < header.n_cols; col++)
        {
            src_properties[row * header.n_cols + col] = *raster_it;
            ++raster_it;
        }
    }

    std::unique_ptr<MeshLib::Mesh> src_mesh(
        MeshLib::RasterToMesh::convert(*raster,
                                       MeshLib::MeshElemType::QUAD,
                                       MeshLib::UseIntensityAs::DATAVECTOR,
                                       property_arg.getValue()));

    // do the interpolation
    MeshLib::Mesh2MeshPropertyInterpolation mesh_interpolation(
        *src_mesh, property_arg.getValue());
    mesh_interpolation.setPropertiesForMesh(dest_mesh.get());

    double scale(1.0);
    if (unit_arg.getValue() == "m/s")
    {
        scale = 1.0;
    }
    else if (unit_arg.getValue() == "mm/a")
    {
        scale = 1e-3 / (365.25 * 86400);
    }
    else if (unit_arg.getValue() == "mm/month")
    {
        scale = 1e-3 * (12.0 / (365.25 * 86400));
    }

    scaleMeshPropertyVector(*dest_mesh, property_arg.getValue(), scale);

    if (!out_mesh_arg.getValue().empty())
    {
        MeshLib::IO::writeMeshToFile(*dest_mesh, out_mesh_arg.getValue());
    }

    return EXIT_SUCCESS;
}
