/*
 * GMSH2OGS.cpp
 *
 *  Created on: Dec 13, 2011
 *      Author: TF
 */

/*
 * mainExtractSurface.cpp
 *
 *  Created on: Jan 26, 2011
 *      Author: TF
 */

// STL
#include <string>

// FileIO
#include "MeshIO/GMSHInterface.h"
#include "MeshIO/OGSMeshIO.h"

// MSH
#include "msh_lib.h" // for FEMRead
#include "msh_mesh.h"

Problem* aproblem = NULL;

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
	std::string file_base_name (tmp);
	if (tmp.find (".msh") != std::string::npos)
		file_base_name = tmp.substr (0, tmp.size() - 4);

	std::vector<MeshLib::CFEMesh*> mesh_vec;
	FEMRead(file_base_name, mesh_vec);
	if (mesh_vec.empty())
	{
		std::cerr << "could not read mesh from file " << tmp << std::endl;
		return -1;
	}
	MeshLib::CFEMesh* mesh (mesh_vec[mesh_vec.size() - 1]);

	// *** create new mesh
	tmp = argv[3];
	if (tmp.find ("--mesh-out") == std::string::npos)
	{
		std::cout << "could not find option --mesh-out" << std::endl;
		return -1;
	}

	if (mesh->GetNodesNumber(false) == 0) {
		mesh->setNumberOfNodesFromNodesVectorSize();
	}

	tmp = argv[4];
	std::cout << "writing mesh to file " << tmp << " ... " << std::flush;
	FileIO::OGSMeshIO mesh_io;
	mesh_io.setMesh(mesh);
	mesh_io.writeToFile (tmp);
	std::cout << "ok" << std::endl;

	delete mesh;

}
