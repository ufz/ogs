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

class FemMesh : public Mesh
{
public:
	FemMesh(const std::string &name, const std::vector<Node*> &nodes, const std::vector<Element*> &elements);
	FemMesh(const Mesh &mesh);
	FemMesh(const FemMesh &mesh);
	FemMesh(const std::string &file_name);
	~FemMesh();


}; /* class */

} /* namespace */

#endif /* FEMMESH_H_ */
