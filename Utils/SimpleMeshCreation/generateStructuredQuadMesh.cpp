/**
 * \file
 * \author Thomas Fischer
 * \date   Sep 21, 2012
 * \brief  Implementation of the generateStructuredQuadMesh tool.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// BaseLib
#include "tclap/CmdLine.h"

// FileIO/Legacy
#include "MeshIO.h"

// Gui/VtkVis
#include "VtkMeshConverter.h"

// MeshLib
#include "Mesh.h"

int main (int argc, char* argv[])
{
	TCLAP::CmdLine cmd("Simple mesh creator for unstructured meshes", ' ', "0.1");
	TCLAP::ValueArg<std::string> mesh_arg("m", "mesh", "the mesh is stored to this file", true, "test.msh", "filename");
	cmd.add( mesh_arg );
	TCLAP::ValueArg<unsigned> width_arg("c","columns","the number of columns of the structured mesh", true, 1000, "number of cols");
	cmd.add( width_arg );
	TCLAP::ValueArg<unsigned> height_arg("r","rows","the number of rows of the structured mesh", true, 1000, "number of rows");
	cmd.add( height_arg );
	TCLAP::ValueArg<double> edge_length_arg("e","edge-length","the size of the edge length", false, 1, "edge length");
	cmd.add( edge_length_arg );
	TCLAP::ValueArg<double> origin_x_arg("x","origin-x","x coordinate of the origin of the mesh", false, 0.0, "x coords");
	cmd.add( origin_x_arg );
	TCLAP::ValueArg<double> origin_y_arg("y","origin-y","y coordinate of the origin of the mesh", false, 0.0, "y coords");
	cmd.add( origin_y_arg );
	TCLAP::ValueArg<double> origin_z_arg("z","origin-z","z coordinate of the origin of the mesh", false, 0.0, "z coords");
	cmd.add( origin_z_arg );

	cmd.parse( argc, argv );

	// thanks to KR for this algorithm
	unsigned height(height_arg.getValue()), width(width_arg.getValue());
	double edge_length(edge_length_arg.getValue());
	const unsigned size (height*width);
	double* values (new double[size]);
	const double origin[3] = {origin_x_arg.getValue() + edge_length/2, origin_y_arg.getValue() + edge_length/2, origin_z_arg.getValue()};
	for (unsigned i=0; i<size; ++i) values[i]=0;
	MeshLib::Mesh* mesh(VtkMeshConverter::convertImgToMesh(values, origin, height, width, edge_length, MeshElemType::QUAD, UseIntensityAs::MATERIAL));

	delete [] values;

	FileIO::Legacy::MeshIO mesh_writer;
	mesh_writer.setMesh(mesh);
	mesh_writer.setPrecision(12);
	mesh_writer.writeToFile(mesh_arg.getValue());

	delete mesh;
}
