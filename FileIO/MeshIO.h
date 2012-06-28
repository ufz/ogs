/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 *
 * \file MeshIO.h
 *
 * Created on 2012-05-08 by Karsten Rink
 */

/**
 * This is currently just test functionality for testing ogs-6 mesh data structures!
 */

#ifndef MESHIO_H_
#define MESHIO_H_

#include "Writer.h"

#include <sstream>
#include <iostream>
#include <vector>

namespace MeshLib
{
	class Mesh;
	class Node;
	class Element;
}

namespace FileIO
{
class MeshIO : public Writer
{
public:
	/// Constructor.
	MeshIO();

	virtual ~MeshIO() {};

	/// Read mesh from file.
	MeshLib::Mesh* loadMeshFromFile(const std::string& fileName);

	/// Set mesh for writing.
	void setMesh(const MeshLib::Mesh*  mesh);

protected:
	/// Write mesh to stream.
	int write(std::ostream &out);

private:
	MeshLib::Element* readElement(const std::string& line, const std::vector<MeshLib::Node*> &nodes);

	double* _edge_length[2];
	const MeshLib::Mesh* _mesh;

};  /* class */

} /* namespace */

#endif /* MESHIO_H_ */

