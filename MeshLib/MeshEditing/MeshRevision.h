/**
 * \file   MeshRevision.h
 * \author Karsten Rink
 * \date   2014-02-14
 * \brief  Definition of the MeshRevision class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MESHREVISION_H_
#define MESHREVISION_H_

#include <array>
#include <vector>
#include <cstddef>

// forward declaration
namespace MeshLib {
	class Mesh;
	class Node;
	class Element;
}

namespace MeshLib {

/**
 * Collapses nodes with a distance smaller min_distance and 
 * reduces elements accordingly.
 */
class MeshRevision {
public:
	/// Constructor
	MeshRevision(Mesh const*const mesh);

	virtual ~MeshRevision() {};

	/**
	 * Collapsed all nodes with distance < eps but ignores elements 
	 * (i.e. elements with collapsed nodes may result)
	 * This is implicitely called when calling simplifyMesh(), so it does not need to be 
	 * called seperately when using simplifyMesh().
	 */
	MeshLib::Mesh* collapseNodes(const std::string &new_mesh_name, double eps);

	/**
	 * Create a new mesh where all nodes with a distance < eps from each other are collapsed.
	 * Elements are adjusted accordingly.
	 * @param eps Minimum distance for nodes not to be collapsed
	 * @param min_elem_dim Minimum dimension of elements to be inserted into new mesh (i.e. min_elem_dim=3 will prevent the new mesh to contain 2D elements)
	 */
	MeshLib::Mesh* simplifyMesh(const std::string &new_mesh_name, double eps, unsigned min_elem_dim = 1);

private:
	/// Designates nodes to be collapsed by setting their ID to the index of the node they will get merged with.
	std::vector<MeshLib::Node*> collapseNodeIndeces(double eps);
	/// Constructs a new node vector for the resulting mesh by removing all nodes whose ID indicates they need to be merged/removed.
	std::vector<MeshLib::Node*> constructNewNodesArray(const std::vector<std::size_t> &id_map);

	/// Calculates the number of unique nodes in an element (i.e. uncollapsed nodes)
	unsigned getNUniqueNodes(MeshLib::Element const*const element) const;

	/// Copies an element without change, using the nodes vector from the result mesh.
	MeshLib::Element* copyElement(MeshLib::Element const*const element, 
		                          const std::vector<MeshLib::Node*> &nodes) const;
	// Revises an element by removing collapsed nodes, using the nodes vector from the result mesh.
	void reduceElement(MeshLib::Element const*const element, 
		               unsigned n_unique_nodes, 
					   const std::vector<MeshLib::Node*> &nodes,
					   std::vector<MeshLib::Element*> &elements,
					   unsigned min_elem_dim) const;

	/// Creates a new line element identical with "line" but using the new nodes vector.
	MeshLib::Element* copyLine(MeshLib::Element const*const line, const std::vector<MeshLib::Node*> &nodes) const;
	/// Creates a new triangle element identical with "tri" but using the new nodes vector.
	MeshLib::Element* copyTri(MeshLib::Element const*const tri, const std::vector<MeshLib::Node*> &nodes) const;
	/// Creates a new quad element identical with "quad" but using the new nodes vector.
	MeshLib::Element* copyQuad(MeshLib::Element const*const quad, const std::vector<MeshLib::Node*> &nodes) const;
	/// Creates a new tetrahedron element identical with "tet" but using the new nodes vector.
	MeshLib::Element* copyTet(MeshLib::Element const*const tet, const std::vector<MeshLib::Node*> &nodes) const;
	/// Creates a new hexahedron element identical with "hex" but using the new nodes vector.
	MeshLib::Element* copyHex(MeshLib::Element const*const hex, const std::vector<MeshLib::Node*> &nodes) const;
	/// Creates a new pyramid element identical with "pyramid" but using the new nodes vector.
	MeshLib::Element* copyPyramid(MeshLib::Element const*const pyramid, const std::vector<MeshLib::Node*> &nodes) const;
	/// Creates a new prism element identical with "prism" but using the new nodes vector.
	MeshLib::Element* copyPrism(MeshLib::Element const*const prism, const std::vector<MeshLib::Node*> &nodes) const;

	/// Creates a line element from the first two unique nodes found in the element (element *should* have exactly two unique nodes!)
	MeshLib::Element* constructLine(MeshLib::Element const*const element, const std::vector<MeshLib::Node*> &nodes) const;
	/// Creates a triangle element from the first three unique nodes found in the element (element *should* have exactly three unique nodes!)
	MeshLib::Element* constructTri(MeshLib::Element const*const element, const std::vector<MeshLib::Node*> &nodes) const;
	/// Creates a quad or a tet, depending if the four nodes being coplanar or not (element *should* have exactly four unique nodes!)
	MeshLib::Element* constructFourNodeElement(MeshLib::Element const*const element, const std::vector<MeshLib::Node*> &nodes, unsigned min_elem_dim = 1) const;

	/** 
	 * Reduces a hexahedron element by removing collapsed nodes and constructing one or more new elements from the remaining nodes.
	 * @return The number of newly created elements
	 */
	unsigned reduceHex(MeshLib::Element const*const hex, 
		           unsigned n_unique_nodes, 
				   const std::vector<MeshLib::Node*> &nodes, 
				   std::vector<MeshLib::Element*> &new_elements,
				   unsigned min_elem_dim) const;
	/// Reduces a pyramid element by removing collapsed nodes and constructing a new elements from the remaining nodes.
	void reducePyramid(MeshLib::Element const*const pyramid, 
		               unsigned n_unique_nodes,
					   const std::vector<MeshLib::Node*> &nodes,
					   std::vector<MeshLib::Element*> &new_elements,
					   unsigned min_elem_dim) const;
	/**
	 * Reduces a prism element by removing collapsed nodes and constructing one or two new elements from the remaining nodes.
	 * @return The number of newly created elements
	 */
	unsigned reducePrism(MeshLib::Element const*const prism,
		             unsigned n_unique_nodes,
					 const std::vector<MeshLib::Node*> &nodes,
					 std::vector<MeshLib::Element*> &new_elements,
					 unsigned min_elem_dim) const;
	
	// In an element with 5 unique nodes, return the node that will be the top of the resulting pyramid
	unsigned findPyramidTopNode(const MeshLib::Element &element, const std::array<unsigned,4> &base_node_ids) const;

	/// Lookup-table for returning the diametral node id of the given node id in a Hex
	unsigned lutHexDiametralNode(unsigned id) const;

	/// Lookup-table for returning four nodes connected to the two nodes (id1, id2) forming an edge in a Hex
	const std::array<unsigned,4> lutHexCuttingQuadNodes(unsigned id1, unsigned id2) const;

	/// When a hex is subdivided into two prisms, this returns the nodes of the hex edge that will serve as the back of one of the prisms.
	const std::pair<unsigned, unsigned> lutHexBackNodes(unsigned i, unsigned j, unsigned k, unsigned l) const;

	/// Lookup-table for returning the third node of bottom or top triangle given the other two
	unsigned lutPrismThirdNode(unsigned id1, unsigned id2) const;

	/// The original mesh used for constructing the class
	Mesh const*const _mesh;
};

}

#endif /* MESHREVISION_H_ */
