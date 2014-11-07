/**
 * \file
 * \author Karsten Rink
 * \date   2013-04-04
 * \brief  Definition of the ElementExtraction
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef ELEMENTEXTRACTION_H
#define ELEMENTEXTRACTION_H

#include <string>
#include <vector>
#include "MeshLib/MeshEnums.h"
#include "MeshLib/Node.h"

namespace MeshLib {

// forward declarations
class Mesh;
class Element;

class ElementExtraction
{
public:
	ElementExtraction(const MeshLib::Mesh &mesh);

	~ElementExtraction();

	/// The error code determined during element extraction
	/// (0 = no errors, 1 = all elements selected, 2 = no elements selected)
	unsigned getErrorCode() { return _error_code; };

	/// Removes all mesh elements marked by search-methods.
	MeshLib::Mesh* removeMeshElements(const std::string &new_mesh_name);

	/// Marks all elements with the given Material ID.
	std::size_t searchByMaterialID(unsigned matID);

	/// Marks all elements of the given element type.
	std::size_t searchByElementType(MeshElemType eleType);

	/// Marks all elements with a volume smaller than eps.
	std::size_t searchByContent(double eps = std::numeric_limits<double>::epsilon());

	/// Marks all elements with at least one node outside the bounding box spanned by x1 and x2;
	std::size_t searchByBoundingBox(const MeshLib::Node &x1, const MeshLib::Node &x2);


private:
	/// Updates the vector of marked elements with values from vec.
	void updateUnion(const std::vector<std::size_t> &vec);

	/// Removes elements from vec_removed in vec_src_elems
	std::vector<MeshLib::Element*> excludeElements(const std::vector<MeshLib::Element*> & vec_src_elems, const std::vector<std::size_t> &vec_removed) const;

	/// The mesh from which elements should be removed.
	const MeshLib::Mesh &_mesh;
	/// The vector of element indices that should be removed.
	std::vector<std::size_t> _marked_elements;
	/// An error code during mesh element extraction for checking the result from outside (0 = no errors).
	unsigned _error_code;
};

} // end namespace MeshLib

#endif //ELEMENTEXTRACTION_H
