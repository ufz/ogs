/**
 * \file
 * \author Thomas Fischer
 * \date   2011-12-13
 * \brief  Implementation of the GMSH2OGS converter.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// STL
#include <string>

// FileIO
#include "Legacy/MeshIO.h"
#include "MeshIO/GMSHInterface.h"

// MeshLib
#include "Mesh.h"

int main (int argc, char* argv[])
{
	if (argc < 5)
	{
		std::cout << "Usage: " << argv[0] <<
		" --mesh-in meshfile --mesh-out meshfile" <<  std::endl;
		return -1;
	}

	// *** read mesh
	std::string tmp (argv[1]);
	if (tmp.find ("--mesh-in") == std::string::npos)
	{
		std::cout << "could not find option --mesh-in" << std::endl;
		return -1;
	}

	tmp = argv[2];
	MeshLib::Mesh* mesh (FileIO::Legacy::MeshIO().loadMeshFromFile(tmp));

	// *** create new mesh
	tmp = argv[3];
	if (tmp.find ("--mesh-out") == std::string::npos)
	{
		std::cout << "could not find option --mesh-out" << std::endl;
		return -1;
	}

	tmp = argv[4];
	std::cout << "writing mesh to file " << tmp << " ... " << std::flush;
	FileIO::Legacy::MeshIO mesh_io;
	mesh_io.setMesh(mesh);
	mesh_io.writeToFile (tmp);
	std::cout << "ok" << std::endl;

	delete mesh;

}
