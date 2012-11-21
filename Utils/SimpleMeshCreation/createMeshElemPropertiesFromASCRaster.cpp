/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 * \file createMeshElemPropertiesFromASCRaster.cpp
 *
 *  Created on  Nov 1, 2012 by Thomas Fischer
 */

// BaseLib
#include "tclap/CmdLine.h"
// ThirdParty/logog
#include "logog/include/logog.hpp"
// BaseLib
#include "StringTools.h"
#include "quicksort.h"
#include "LogogSimpleFormatter.h"

// FileIO/Legacy
#include "MeshIO.h"
#include "readMeshFromFile.h"
#include "FileTools.h"

// GeoLib
#include "Raster.h"

// Gui/VtkVis
#include "VtkMeshConverter.h"

// MeshLib
#include "Mesh.h"
#include "Elements/Element.h"
#include "MshEnums.h"
#include "Mesh2MeshPropertyInterpolation.h"
#include "ConvertRasterToMesh.h"

int main (int argc, char* argv[])
{
	LOGOG_INITIALIZE();
	logog::Cout* logog_cout (new logog::Cout);
	BaseLib::LogogSimpleFormatter *custom_format (new BaseLib::LogogSimpleFormatter);
	logog_cout->SetFormatter(*custom_format);

	TCLAP::CmdLine cmd("Generates properties for mesh elements of an input mesh deploying a ASC raster file", ' ', "0.1");

	TCLAP::ValueArg<std::string> out_mesh_arg("o", "out-mesh", "the mesh is stored to a file of this name", false, "", "filename for mesh output");
	cmd.add( out_mesh_arg );

	TCLAP::ValueArg<bool> refinement_raster_output_arg("","output-refined-raster", "write refined raster to a new ASC file", false, false, "0");
	cmd.add( refinement_raster_output_arg );

	TCLAP::ValueArg<unsigned> refinement_arg("r","refine","refinement factor that raises the resolution of the raster data", false, 1, "factor (default = 1)");
	cmd.add( refinement_arg );

	TCLAP::ValueArg<std::string> raster_arg("","raster-file","the name of the ASC raster file", true, "", "file name");
	cmd.add( raster_arg );

	TCLAP::ValueArg<std::string> mesh_arg("m", "mesh", "the mesh is read from this file", true, "test.msh", "filename");
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
			std::string new_raster_fname (BaseLib::dropFileExtension(raster_arg.getValue()));
			new_raster_fname += "-" + BaseLib::number2str(raster->getNRows()) + "x" + BaseLib::number2str(raster->getNCols()) + ".asc";
			std::ofstream out(new_raster_fname);
			raster->writeRasterAsASC(out);
			out.close();
		}
	}

	// put raster data in a std::vector
	GeoLib::Raster::const_iterator raster_it(raster->begin());
	unsigned n_cols(raster->getNCols()), n_rows(raster->getNRows());
	std::vector<double> src_properties(n_cols*n_rows);
	for (unsigned row(0); row<n_rows; row++) {
		for (unsigned col(0); col<n_cols; col++) {
			src_properties[row*n_cols+col] = *raster_it;
			++raster_it;
		}
	}

	{
		double src_mean_value(src_properties[0]);
		for (size_t k(1); k<n_cols*n_rows; k++) {
			src_mean_value += src_properties[k];
		}
		src_mean_value /= n_cols*n_rows;
		std::cout << "mean value of source: " << src_mean_value << std::endl;

		double src_varianz(MathLib::fastpow(src_properties[0] - src_mean_value, 2));
		for (size_t k(1); k<n_cols*n_rows; k++) {
			src_varianz += MathLib::fastpow(src_properties[k] - src_mean_value, 2);
		}
		src_varianz /= n_cols*n_rows;
		std::cout << "variance of source: " << src_varianz << std::endl;
	}

//	double spacing(raster->getRasterPixelDistance());
//	double *raster_with_alpha(new double[raster->getNRows() * raster->getNCols()]);
//	raster_it = raster->begin();
//	for (std::size_t k(0); k<raster->getNRows() * raster->getNCols(); k++) {
//		raster_with_alpha[k] = *raster_it;
//		++raster_it;
//	}

//	double origin[3] = {raster->getOrigin()[0] + spacing/2.0, raster->getOrigin()[1] + spacing/2.0, raster->getOrigin()[2]};
//	MeshLib::Mesh* src_mesh (VtkMeshConverter::convertImgToMesh(raster_with_alpha, origin, raster->getNCols(), raster->getNRows(),
//					spacing, MshElemType::QUAD, UseIntensityAs::MATERIAL));

	MeshLib::Mesh* src_mesh(MeshLib::ConvertRasterToMesh(*raster, MshElemType::QUAD,
					MeshLib::UseIntensityAs::MATERIAL).execute());

//	delete [] raster_with_alpha;

	std::vector<size_t> src_perm(n_cols*n_rows);
	for (size_t k(0); k<n_cols*n_rows; k++) src_perm[k] = k;
	BaseLib::Quicksort<double>(src_properties, 0, n_cols*n_rows, src_perm);

	// compress the property data structure
	const size_t mat_map_size(src_properties.size());
	std::vector<size_t> mat_map(mat_map_size);
	mat_map[0] = 0;
	size_t n_mat(1);
	for (size_t k(1); k<mat_map_size; ++k) {
		if (std::fabs(src_properties[k] - src_properties[k-1]) > std::numeric_limits<double>::epsilon()) {
			mat_map[k] = mat_map[k-1]+1;
			n_mat++;
		} else
			mat_map[k] = mat_map[k-1];
	}
	std::vector<double> compressed_src_properties(n_mat);
	compressed_src_properties[0] = src_properties[0];
	for (size_t k(1), id(1); k<mat_map_size; ++k) {
		if (std::fabs(src_properties[k] - src_properties[k-1]) > std::numeric_limits<double>::epsilon()) {
			compressed_src_properties[id] = src_properties[k];
			id++;
		}
	}
	compressed_src_properties[n_mat-1] = src_properties[mat_map_size-1];

	// reset materials in source mesh
	const size_t n_mesh_elements(src_mesh->getNElements());
	for (size_t k(0); k<n_mesh_elements; k++) {
		const_cast<MeshLib::Element*>(src_mesh->getElement(src_perm[k]))->setValue(mat_map[k]);
	}

	// do the interpolation
	MeshLib::Mesh2MeshPropertyInterpolation mesh_interpolation(src_mesh, &compressed_src_properties);
	std::vector<double> dest_properties(dest_mesh->getNElements());
	mesh_interpolation.setPropertiesForMesh(const_cast<MeshLib::Mesh*>(dest_mesh), dest_properties);

	const size_t n_dest_mesh_elements(dest_mesh->getNElements());

	{ // write property file
		std::string property_fname(BaseLib::dropFileExtension(raster_arg.getValue())+".hd");
		std::ofstream property_out(property_fname.c_str());
		if (! property_out) {
			std::cerr << "could not open file " << property_fname << "PropertyMapping" << std::endl;
			return -1;
		}

		for (size_t k(0); k<n_dest_mesh_elements; k++) {
			property_out << k << " " << dest_properties[k] << std::endl;
		}
		property_out.close();
	}

	{
		double mu(dest_properties[0]);
		for (size_t k(1); k<n_dest_mesh_elements; k++) {
			mu += dest_properties[k];
		}
		mu /= n_dest_mesh_elements;
		std::cout << "mean value of destination: " << mu << std::endl;

		double sigma_q(MathLib::fastpow(dest_properties[0] - mu, 2));
		for (size_t k(1); k<n_dest_mesh_elements; k++) {
			sigma_q += MathLib::fastpow(dest_properties[k] - mu, 2);
		}
		sigma_q /= n_dest_mesh_elements;
		std::cout << "variance of destination: " << sigma_q << std::endl;
	}

	if (! out_mesh_arg.getValue().empty()) {
		std::vector<size_t> dest_perm(n_dest_mesh_elements);
		for (size_t k(0); k<n_dest_mesh_elements; k++) dest_perm[k] = k;
		BaseLib::Quicksort<double>(dest_properties, 0, n_dest_mesh_elements, dest_perm);

		// reset materials in destination mesh
		for (size_t k(0); k<n_dest_mesh_elements; k++) {
			const_cast<MeshLib::Element*>(dest_mesh->getElement(dest_perm[k]))->setValue(k);
		}

		FileIO::MeshIO mesh_writer;
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
