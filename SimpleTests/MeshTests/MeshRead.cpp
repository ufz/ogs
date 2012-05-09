/*
 * MeshRead.cpp
 *
 *  Created on: 2012/05/09
 *      Author: KR
 */

#include "Mesh.h"
#include "MeshIO.h"

int main(int argc, char *argv[])
{
	std::string file_name("c:/Project/Data/Ammer/Ammer-Homogen100m-Final.msh");
	//std::string file_name("c:/Project/PlyTestMesh.msh");
	FileIO::MeshIO mesh_io;
	MeshLib::Mesh* mesh = mesh_io.loadMeshFromFile(file_name);

	delete mesh;
}

