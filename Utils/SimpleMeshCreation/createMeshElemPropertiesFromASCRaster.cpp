/**
 * \file
 * \author Thomas Fischer
 * \date   Nov 1, 2012
 * \brief  Implementation of the createMeshElemPropertiesFromASCRaster tool.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// stl
#include <numeric>

// BaseLib
#include "tclap/CmdLine.h"

// ThirdParty/logog
#include "logog/include/logog.hpp"

// BaseLib
#include "LogogSimpleFormatter.h"
#include "quicksort.h"

// FileIO/Legacy
#include "FileTools.h"
#include "MeshIO.h"
#include "readMeshFromFile.h"
#include "AsciiRasterInterface.h"

// GeoLib
#include "Raster.h"

// Gui/VtkVis
#include "VtkMeshConverter.h"

// MathLib
#include "MathTools.h"

// MeshLib
#include "MeshGenerators/ConvertRasterToMesh.h"
#include "Elements/Element.h"
#include "Mesh.h"
#include "MeshEditing/Mesh2MeshPropertyInterpolation.h"
#include "MeshEnums.h"

int main (int argc, char* argv[])
{
	LOGOG_INITIALIZE();
	logog::Cout* logog_cout (new logog::Cout);
	BaseLib::LogogSimpleFormatter *custom_format (new BaseLib::LogogSimpleFormatter);
	logog_cout->SetFormatter(*custom_format);

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
	MeshLib::Mesh* dest_mesh(FileIO::readMeshFromFile(mesh_arg.getValue()));

	// read raster and if required manipulate it
	GeoLib::Raster* raster(GeoLib::Raster::getRasterFromASCFile(raster_arg.getValue()));
	if (refinement_arg.getValue() > 1) {
		raster->refineRaster(refinement_arg.getValue());
		if (refinement_raster_output_arg.getValue()) {
			// write new asc file
			std::string new_raster_fname (BaseLib::dropFileExtension(
			                                      raster_arg.getValue()));
			new_raster_fname += "-" + std::to_string(raster->getNRows()) + "x" +
			                    std::to_string(raster->getNCols()) + ".asc";
			FileIO::AsciiRasterInterface::writeRasterAsASC(raster, new_raster_fname);
		}
	}

	// put raster data in a std::vector
	GeoLib::Raster::const_iterator raster_it(raster->begin());
	std::size_t n_cols(raster->getNCols()), n_rows(raster->getNRows());
	std::size_t size(n_cols * n_rows);
	std::vector<double> src_properties(size);
	for (unsigned row(0); row<n_rows; row++) {
		for (unsigned col(0); col<n_cols; col++) {
			src_properties[row * n_cols + col] = *raster_it;
			++raster_it;
		}
	}

	{
		const double mu(std::accumulate(src_properties.begin(), src_properties.end(), 0.0) / size);
		INFO("Mean value of source: %f.", mu);

		double src_variance(MathLib::fastpow(src_properties[0] - mu, 2));
		for (std::size_t k(1); k<size; k++) {
			src_variance += MathLib::fastpow(src_properties[k] - mu, 2);
		}
		src_variance /= size;
		INFO("Variance of source: %f.", src_variance);
	}

	MeshLib::Mesh* src_mesh(MeshLib::ConvertRasterToMesh(*raster, MeshElemType::QUAD,
					MeshLib::UseIntensityAs::MATERIAL).execute());

	std::vector<std::size_t> src_perm(size);
	std::iota(src_perm.begin(), src_perm.end(), 0);
	BaseLib::Quicksort<double>(src_properties, 0, size, src_perm);

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
	const std::size_t n_mesh_elements(src_mesh->getNElements());
	for (std::size_t k(0); k<n_mesh_elements; k++) {
		const_cast<MeshLib::Element*>(src_mesh->getElement(src_perm[k]))->setValue(mat_map[k]);
	}

	// do the interpolation
	MeshLib::Mesh2MeshPropertyInterpolation mesh_interpolation(src_mesh,
	                                                           &compressed_src_properties);
	std::vector<double> dest_properties(dest_mesh->getNElements());
	mesh_interpolation.setPropertiesForMesh(const_cast<MeshLib::Mesh*>(dest_mesh),
	                                        dest_properties);

	const std::size_t n_dest_mesh_elements(dest_mesh->getNElements());

	{ // write property file
		std::string property_fname(mapping_arg.getValue());
		std::ofstream property_out(property_fname.c_str());
		if (!property_out)
		{
			ERR("Could not open file %s for writing the mapping.", property_fname.c_str());
			return -1;
		}

		for (std::size_t k(0); k < n_dest_mesh_elements; k++)
			property_out << k << " " << dest_properties[k] << "\n";
		property_out.close();
	}

	{
		const double mu(std::accumulate(dest_properties.begin(), dest_properties.end(), 0.0) / n_dest_mesh_elements);
		INFO("Mean value of destination: %f.", mu);

		double sigma_q(MathLib::fastpow(dest_properties[0] - mu, 2));
		for (std::size_t k(1); k < n_dest_mesh_elements; k++)
			sigma_q += MathLib::fastpow(dest_properties[k] - mu, 2);
		sigma_q /= n_dest_mesh_elements;
		INFO("Variance of destination: %f.", sigma_q);
	}

	if (! out_mesh_arg.getValue().empty()) {
		std::vector<std::size_t> dest_perm(n_dest_mesh_elements);
		std::iota(dest_perm.begin(), dest_perm.end(), 0);
		BaseLib::Quicksort<double>(dest_properties, 0, n_dest_mesh_elements, dest_perm);

		// reset materials in destination mesh
		for (std::size_t k(0); k<n_dest_mesh_elements; k++) {
			const_cast<MeshLib::Element*>(dest_mesh->getElement(dest_perm[k]))->setValue(k);
		}

		FileIO::Legacy::MeshIO mesh_writer;
		mesh_writer.setPrecision(12);
		mesh_writer.setMesh(dest_mesh);
		mesh_writer.writeToFile(out_mesh_arg.getValue());
	}

	delete raster;
	delete src_mesh;
	delete dest_mesh;

	delete custom_format;
	delete logog_cout;
	LOGOG_SHUTDOWN();

	return 0;
}
