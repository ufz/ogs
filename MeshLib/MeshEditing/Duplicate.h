/**
 * \file   Duplicate.h
 * \author Karsten Rink
 * \date   2014-03-25
 * \brief  Definition of Duplicate functions
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef DUPLICATE_H
#define DUPLICATE_H

#include <vector>

namespace MeshLib 
{

	class Mesh;
	class Node;
	class Element;

namespace Duplicate 
{

	/// Creates a deep copy of a Node vector
	std::vector<MeshLib::Node*> NodeVector(const std::vector<MeshLib::Node*> &nodes);

	/// Creates a deep copy of an element vector using the given Node vector
	std::vector<MeshLib::Element*> NodeVector(const std::vector<MeshLib::Element*> &elements, const std::vector<MeshLib::Node*> &nodes);

	/// Copies an element without change, using the nodes vector from the result mesh.
	MeshLib::Element* copyElement(MeshLib::Element const*const element, 
		                          const std::vector<MeshLib::Node*> &nodes);

	/// Creates a new line element identical with "line" but using the new nodes vector.
	MeshLib::Element* copyLine(MeshLib::Element const*const line, const std::vector<MeshLib::Node*> &nodes);
	/// Creates a new triangle element identical with "tri" but using the new nodes vector.
	MeshLib::Element* copyTri(MeshLib::Element const*const tri, const std::vector<MeshLib::Node*> &nodes);
	/// Creates a new quad element identical with "quad" but using the new nodes vector.
	MeshLib::Element* copyQuad(MeshLib::Element const*const quad, const std::vector<MeshLib::Node*> &nodes);
	/// Creates a new tetrahedron element identical with "tet" but using the new nodes vector.
	MeshLib::Element* copyTet(MeshLib::Element const*const tet, const std::vector<MeshLib::Node*> &nodes);
	/// Creates a new hexahedron element identical with "hex" but using the new nodes vector.
	MeshLib::Element* copyHex(MeshLib::Element const*const hex, const std::vector<MeshLib::Node*> &nodes);
	/// Creates a new pyramid element identical with "pyramid" but using the new nodes vector.
	MeshLib::Element* copyPyramid(MeshLib::Element const*const pyramid, const std::vector<MeshLib::Node*> &nodes);
	/// Creates a new prism element identical with "prism" but using the new nodes vector.
	MeshLib::Element* copyPrism(MeshLib::Element const*const prism, const std::vector<MeshLib::Node*> &nodes);

} // end namespace Duplicate

} // end namespace MeshLib

#endif //DUPLICATE_H
