/*
 * MeshRead.cpp
 *
 *  Created on: 2012/05/09
 *      Author: KR
 */

// BaseLib
#include "MemWatch.h"
#include "RunTime.h"

// MeshLib
#include "Node.h"
#include "Elements/Element.h"
#include "Mesh.h"
#include "MeshIO.h"

int main(int argc, char *argv[])
{
	//std::string file_name("/mnt/visdata/tom/data/TestMeshes/Mesh-dx1.00-Layered20.msh");
	std::string file_name("c:/Project/PlyTestMesh.msh");
	std::cout << "sizeof(double): " << sizeof (double) << std::endl;
	std::cout << "sizeof(GeoLib::Point): " << sizeof (GeoLib::Point) << std::endl;
	std::cout << "sizeof(GeoLib::PointWithID): " << sizeof (GeoLib::PointWithID) << std::endl;
	std::cout << "sizeof(Node): " << sizeof (MeshLib::Node) << std::endl;
	std::cout << "sizeof(Element): " << sizeof (MeshLib::Element) << std::endl;
	FileIO::MeshIO mesh_io;
#ifndef WIN32
	BaseLib::MemWatch mem_watch;
	unsigned long mem_without_mesh (mem_watch.getVirtMemUsage());
	BaseLib::RunTime run_time;
	run_time.start();
#endif
	MeshLib::Mesh* mesh = mesh_io.loadMeshFromFile(file_name);
#ifndef WIN32
	unsigned long mem_with_mesh (mem_watch.getVirtMemUsage());
	std::cout << "mem for mesh: " << (mem_with_mesh - mem_without_mesh)/(1024*1024) << " MB" << std::endl;
	run_time.stop();
	std::cout << "time for reading: " << run_time.elapsed() << " s" << std::endl;
#endif
	delete mesh;
}

