/**
 * \file
 * \author Karsten Rink
 * \date   2012-05-08
 * \brief  Definition of the MeshIO class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */


#ifndef MESHIO_H_
#define MESHIO_H_

#include "Writer.h"
#include "MeshEnums.h"

#include <sstream>
#include <iostream>
#include <string>
#include <vector>

namespace MeshLib
{
	class Mesh;
	class Node;
	class Element;
}


namespace FileIO
{
namespace Legacy {

/// Interface for handling mesh files from OGS-5 and below. (*.msh files)
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
	void writeElements(std::vector<MeshLib::Element*> const& ele_vec, std::ostream &out) const;
	MeshLib::Element* readElement(const std::string& line, const std::vector<MeshLib::Node*> &nodes);
	const std::string ElemType2StringOutput(const MeshElemType t) const;

	double* _edge_length[2];
	const MeshLib::Mesh* _mesh;

};  /* class */

}
} /* namespace FileIO */

#endif /* MESHIO_H_ */

