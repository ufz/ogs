/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef ELEMENTSEARCH_H_
#define ELEMENTSEARCH_H_

#include <limits>
#include <vector>

#include "GeoLib/AABB.h"
#include "MeshLib/MeshEnums.h"

namespace MeshLib {

// forward declarations
class Mesh;
class Element;

/// Element search class
class ElementSearch final
{
public:
	explicit ElementSearch(const MeshLib::Mesh &mesh);

	/// return marked elements
	const std::vector<std::size_t>& getSearchedElementIDs() const { return _marked_elements; }

	/// Marks all elements with the given Material ID.
	std::size_t searchByMaterialID(int const matID);

	/// Marks all elements of the given element type.
	std::size_t searchByElementType(MeshElemType eleType);

	/// Marks all elements with a volume smaller than eps.
	std::size_t searchByContent(double eps = std::numeric_limits<double>::epsilon());

	/// Marks all elements with at least one node outside the bounding box spanned by x1 and x2;
	std::size_t searchByBoundingBox(GeoLib::AABB const& aabb);

	/// Marks all elements connecting to any of the given nodes
	std::size_t searchByNodeIDs(const std::vector<std::size_t> &node_ids);

private:
	/// Updates the vector of marked elements with values from vec.
	void updateUnion(const std::vector<std::size_t> &vec);

	/// The mesh from which elements should be removed.
	const MeshLib::Mesh &_mesh;
	/// The vector of element indices that should be removed.
	std::vector<std::size_t> _marked_elements;
};

} // end namespace MeshLib

#endif //ELEMENTEXTRACTION_H
