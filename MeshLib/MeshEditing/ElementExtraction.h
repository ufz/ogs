/**
 * \file
 * \author Karsten Rink
 * \date   2013-04-04
 * \brief  Definition of the ElementExtraction
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef ELEMENTEXTRACTION_H
#define ELEMENTEXTRACTION_H

#include <string>
#include <vector>
#include "MeshEnums.h"
#include "Node.h"

namespace MeshLib {

// forward declarations
class Mesh;
class Element;

class ElementExtraction
{
public:
	ElementExtraction(const MeshLib::Mesh &mesh);

	~ElementExtraction();

	/// Removes all mesh elements marked by search-methods.
	MeshLib::Mesh* removeMeshElements(const std::string &new_mesh_name) const;
	
	/// Marks all elements with the given Material ID.
	void searchByMaterialID(unsigned matID);

	/// Marks all elements of the given element type.
	void searchByElementType(MeshElemType eleType);

	/// Marks all elements with a volume smaller than std::numeric_limits<double>::epsilon().
	void searchByZeroContent();

	/// Marks all elements with at least one node outside the bounding box spanned by x1 and x2;
	void searchByBoundingBox(const MeshLib::Node &x1, const MeshLib::Node &x2);
	

private:
	/// Updates the vector of marked elements with values from vec.
	void updateUnion(const std::vector<std::size_t> &vec);

	/// Removes elements from vec_removed in vec_src_elems
	std::vector<MeshLib::Element*> excludeElements(const std::vector<MeshLib::Element*> & vec_src_elems, const std::vector<std::size_t> &vec_removed) const;

	/// Copies nodes and elements of the original mesh for constructing the new mesh
	void copyNodesElements(const std::vector<MeshLib::Node*> &src_nodes, const std::vector<MeshLib::Element*> &src_elems,
						   std::vector<MeshLib::Node*> &dst_nodes, std::vector<MeshLib::Element*> &dst_elems) const;

	const MeshLib::Mesh &_mesh;
	std::vector<std::size_t> _marked_elements;
};

} // end namespace MeshLib

#endif //ELEMENTEXTRACTION_H
