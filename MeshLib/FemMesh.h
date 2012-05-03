/**
 * FemMesh.h
 *
 *      Date: 2012/05/02
 *      Author: KR
 */

#ifndef FEMMESH_H_
#define FEMMESH_H_

#include "Mesh.h"

namespace MeshLib {

/**
 * A finite element mesh.
 */
class FemMesh : public Mesh
{
public:
	/// Constructor using a mesh name and an array of nodes and elements
	FemMesh(const std::string &name, const std::vector<Node*> &nodes, const std::vector<Element*> &elements);

	/// Constructor using a basic mesh
	FemMesh(const Mesh &mesh);

	/// Copy constructor
	FemMesh(const FemMesh &mesh);

	/// Constructor for reading a mesh from a file
	FemMesh(const std::string &file_name);

	/// Destructor
	~FemMesh();


}; /* class */

} /* namespace */

#endif /* FEMMESH_H_ */

