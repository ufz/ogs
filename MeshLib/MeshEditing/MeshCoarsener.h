/**
 * \file
 * \author Thomas Fischer
 * \date   Aug 3, 2012
 * \brief  Definition of the MeshCoarsener class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MESHCOARSENER_H_
#define MESHCOARSENER_H_

#include <vector>
#include <cstddef>

// forward declaration
namespace MeshLib {
	class Mesh;
	class Node;
	class Element;
	class Tri;
	class Quad;
	class Tet;
	class Hex;
	class Pyramid;
	class Prism;
}

namespace MeshLib {

/**
 * The class MeshCoarsener merges mesh nodes that have a smaller
 * distance than a (user) given minimal distance.
 */
class MeshCoarsener {
public:
	/**
	 * Constructor of class MeshCoarsener that takes the mesh object that should be coarsened.
	 * @param mesh the mesh object
	 */
	MeshCoarsener(Mesh const*const mesh);

	/**
	 * destructor.
	 */
	virtual ~MeshCoarsener();

	/**
	 * create new mesh and apply the coarsening process to the mesh
	 * @param min_distance
	 */
	Mesh* operator() (double min_distance);

	Mesh* simplifyMesh();

private:

	std::vector<MeshLib::Node*> collapseNodes(std::vector<unsigned> &idx_map, double eps);

	unsigned getNUniqueNodes(const MeshLib::Element* element);

	MeshLib::Element* reduceTri(MeshLib::Element* tri, unsigned unique_nodes);
	MeshLib::Element* reduceQuad(MeshLib::Element* quad, unsigned unique_nodes);
	MeshLib::Element* reduceTet(MeshLib::Element* tet, unsigned unique_nodes);
	MeshLib::Element* reduceHex(MeshLib::Element* hex, unsigned unique_nodes);
	MeshLib::Element* reducePyramid(MeshLib::Element* pyramid, unsigned unique_nodes);
	MeshLib::Element* reducePrism(MeshLib::Element* prism, unsigned unique_nodes);

	Mesh const*const _orig_mesh;
};

}

#endif /* MESHCOARSENER_H_ */
