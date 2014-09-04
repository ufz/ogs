/**
 * @copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */
#ifndef BOUNDARYELEMENTSSEARCHER_H_
#define BOUNDARYELEMENTSSEARCHER_H_

#include <vector>

// GeoLib
#include "Point.h"
#include "Polyline.h"

// MeshLib
#include "Mesh.h"
#include "Node.h"

// MeshGeoToolsLib
#include "MeshNodeSearcher.h"

namespace MeshGeoToolsLib
{

/**
 * This class searches and creates boundary elements located on a given geometric object.
 * Boundary elements will be created from edges or faces of existing domain elements.
 */
class BoundaryElementsSearcher
{
public:
	/**
	 * Constructor
	 * @param mesh             a mesh object
	 * @param mshNodeSearcher  a MeshNodeSearcher object which is internally used to search mesh nodes
	 */
	BoundaryElementsSearcher(MeshLib::Mesh const& mesh, MeshNodeSearcher &mshNodeSearcher);

	/// destructor
	virtual ~BoundaryElementsSearcher() {}

	/**
	 * generate boundary elements on the given geometric object (point, polyline, surface).
	 *
	 * @param geoObj a GeoLib::GeoObject where the nearest mesh node is searched for
	 * @return a vector of boundary element objects
	 */
	std::vector<MeshLib::Element*> getBoundaryElements(GeoLib::GeoObject const& geoObj);

	/**
	 * generate boundary elements on the given polyline.
	 * @param ply the GeoLib::Polyline the nearest mesh nodes are searched for
	 * @return a vector of boundary element objects
	 */
	std::vector<MeshLib::Element*> getBoundaryElementsAlongPolyline(GeoLib::Polyline const& ply);

private:
	MeshLib::Mesh const& _mesh;
	MeshNodeSearcher &_mshNodeSearcher;
};

} // end namespace MeshGeoTools

#endif /* BOUNDARYELEMENTSSEARCHER_H_ */
