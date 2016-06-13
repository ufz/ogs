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
#include <numeric>

#include <tclap/CmdLine.h>

#include "Applications/ApplicationsLib/LogogSetup.h"

#include "BaseLib/quicksort.h"
#include "BaseLib/FileTools.h"

#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/IO/writeMeshToFile.h"
#include "GeoLib/IO/AsciiRasterInterface.h"

#include "GeoLib/Raster.h"

#include "MathLib/MathTools.h"

#include "MeshLib/MeshGenerators/RasterToMesh.h"
#include "MeshLib/MeshGenerators/VtkMeshConverter.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshEditing/Mesh2MeshPropertyInterpolation.h"
#include "MeshLib/MeshEnums.h"

// From wikipedia:
// http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance. The
// original is citing D.E. Knuth. TAOCP, vol 2.
template <typename InputIterator>
auto computeMeanAndVariance(InputIterator first, InputIterator last) ->
    std::pair<typename InputIterator::value_type, typename InputIterator::value_type>
{
    using T = typename InputIterator::value_type;
    std::size_t n = 0;
    auto mu = T{0};
    auto M2 = T{0};

    while (first != last)
    {
        T const x = *first++;
        n++;
        auto delta = x - mu;
        mu += delta/n;
        M2 += delta * (x - mu);
    }

    if (n < 2)
        return std::make_pair(mu, T{0});

    return std::make_pair(mu, M2/(n - 1));
}

int main (int argc, char* argv[])
{
    ApplicationsLib::LogogSetup logo_setup;

    TCLAP::CmdLine cmd(
            "Generates properties for mesh elements of an input mesh deploying a ASC raster file",
            ' ',
            "0.1");

    TCLAP::ValueArg<std::string> out_mesh_arg("o",
                                              "out-mesh",
                                              "the mesh is stored to a file of this name",
                                              false,
                                              "",
                                              "filename for mesh output");
    cmd.add( out_mesh_arg );

    TCLAP::ValueArg<bool> refinement_raster_output_arg("",
                                                       "output-refined-raster",
                                                       "write refined raster to a new ASC file",
                                                       false,
                                                       false,
                                                       "0");
    cmd.add( refinement_raster_output_arg );

    TCLAP::ValueArg<unsigned> refinement_arg(
            "r",
            "refine",
            "refinement factor that raises the resolution of the raster data",
            false,
            1,
            "factor (default = 1)");
    cmd.add( refinement_arg );

    TCLAP::ValueArg<std::string> mapping_arg("",
                                             "mapping-name",
                                             "file name of mapping",
                                             true,
                                             "",
                                             "file name");
    cmd.add( mapping_arg );

    TCLAP::ValueArg<std::string> raster_arg("",
                                            "raster-file",
                                            "the name of the ASC raster file",
                                            true,
                                            "",
                                            "file name");
    cmd.add( raster_arg );

    TCLAP::ValueArg<std::string> mesh_arg("m",
                                          "mesh",
                                          "the mesh is read from this file",
                                          true,
                                          "test.msh",
                                          "file name");
    cmd.add( mesh_arg );

    cmd.parse( argc, argv );

    // read mesh
    std::unique_ptr<MeshLib::Mesh> dest_mesh(MeshLib::IO::readMeshFromFile(mesh_arg.getValue()));

    // read raster and if required manipulate it
    auto raster = std::unique_ptr<GeoLib::Raster>(
        GeoLib::IO::AsciiRasterInterface::getRasterFromASCFile(raster_arg.getValue()));
    GeoLib::RasterHeader header (raster->getHeader());
    if (refinement_arg.getValue() > 1) {
        raster->refineRaster(refinement_arg.getValue());
        if (refinement_raster_output_arg.getValue()) {
            // write new asc file
            std::string new_raster_fname (BaseLib::dropFileExtension(
                                                  raster_arg.getValue()));
            new_raster_fname += "-" + std::to_string(header.n_rows) + "x" +
                                std::to_string(header.n_cols) + ".asc";
            GeoLib::IO::AsciiRasterInterface::writeRasterAsASC(*raster, new_raster_fname);
        }
    }

    // put raster data in a std::vector
    GeoLib::Raster::const_iterator raster_it(raster->begin());
    std::size_t size(header.n_cols * header.n_rows);
    std::vector<double> src_properties(size);
    for (unsigned row(0); row<header.n_rows; row++) {
        for (unsigned col(0); col<header.n_cols; col++) {
            src_properties[row * header.n_cols + col] = *raster_it;
            ++raster_it;
        }
    }

    {
        double mu, var;
        std::tie(mu, var) = computeMeanAndVariance(src_properties.begin(), src_properties.end());
        INFO("Mean value of source: %f.", mu);
        INFO("Variance of source: %f.", var);
    }

    std::unique_ptr<MeshLib::Mesh> src_mesh(MeshLib::RasterToMesh::convert(
        *raster, MeshLib::MeshElemType::QUAD,MeshLib::UseIntensityAs::DATAVECTOR));

    std::vector<std::size_t> src_perm(size);
    std::iota(src_perm.begin(), src_perm.end(), 0);
    BaseLib::quicksort<double>(src_properties, 0, size, src_perm);

    // compress the property data structure
    const std::size_t mat_map_size(src_properties.size());
    std::vector<std::size_t> mat_map(mat_map_size);
    mat_map[0] = 0;
    std::size_t n_mat(1);
    for (std::size_t k(1); k<mat_map_size; ++k) {
        if (std::fabs(src_properties[k] - src_properties[k-1]) > std::numeric_limits<double>::epsilon()) {
            mat_map[k] = mat_map[k - 1] + 1;
            n_mat++;
        } else
            mat_map[k] = mat_map[k - 1];
    }
    std::vector<double> compressed_src_properties(n_mat);
    compressed_src_properties[0] = src_properties[0];
    for (std::size_t k(1), id(1); k<mat_map_size; ++k) {
        if (std::fabs(src_properties[k] - src_properties[k-1]) > std::numeric_limits<double>::epsilon()) {
            compressed_src_properties[id] = src_properties[k];
            id++;
        }
    }
    compressed_src_properties[n_mat - 1] = src_properties[mat_map_size - 1];

    // reset materials in source mesh
    const std::size_t n_mesh_elements(src_mesh->getNumberOfElements());
    auto materialIds = src_mesh->getProperties().getPropertyVector<int>("MaterialIDs");
    if (!materialIds)
    {
        materialIds = src_mesh->getProperties().createNewPropertyVector<int>
            ("MaterialIDs", MeshLib::MeshItemType::Cell, 1);
        materialIds->resize(n_mesh_elements, 0);    // default material id
    }
    for (std::size_t k(0); k<n_mesh_elements; k++)
        (*materialIds)[src_mesh->getElement(src_perm[k])->getID()] = mat_map[k];

    // do the interpolation
    MeshLib::Mesh2MeshPropertyInterpolation mesh_interpolation(src_mesh.get(),
                                                               &compressed_src_properties);
    std::vector<double> dest_properties(dest_mesh->getNumberOfElements());
    mesh_interpolation.setPropertiesForMesh(dest_mesh.get(),
                                            dest_properties);

    const std::size_t n_dest_mesh_elements(dest_mesh->getNumberOfElements());

    { // write property file
        std::string property_fname(mapping_arg.getValue());
        std::ofstream property_out(property_fname.c_str());
        if (!property_out)
        {
            ERR("Could not open file %s for writing the mapping.", property_fname.c_str());
            return EXIT_FAILURE;
        }

        for (std::size_t k(0); k < n_dest_mesh_elements; k++)
            property_out << k << " " << dest_properties[k] << "\n";
        property_out.close();
    }

    {
        double mu, var;
        std::tie(mu, var) = computeMeanAndVariance(dest_properties.begin(), dest_properties.end());
        INFO("Mean value of destination: %f.", mu);
        INFO("Variance of destination: %f.", var);
    }

    if (! out_mesh_arg.getValue().empty()) {
        std::vector<std::size_t> dest_perm(n_dest_mesh_elements);
        std::iota(dest_perm.begin(), dest_perm.end(), 0);
        BaseLib::quicksort<double>(dest_properties, 0, n_dest_mesh_elements, dest_perm);

        // reset materials in destination mesh
        materialIds = dest_mesh->getProperties().getPropertyVector<int>("MaterialIDs");
        for (std::size_t k(0); k<n_dest_mesh_elements; k++) {
            (*materialIds)[dest_mesh->getElement(dest_perm[k])->getID()] = k;
        }

        MeshLib::IO::writeMeshToFile(*dest_mesh, out_mesh_arg.getValue());
    }

    return EXIT_SUCCESS;
}
