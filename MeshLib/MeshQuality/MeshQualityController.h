/**
 * \file   MeshQualityController.h
 * \author Karsten Rink
 * \date   2013-04-04
 * \brief  Definition of the MeshQualityController class
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MESHQUALITYCONTROLLER_H
#define MESHQUALITYCONTROLLER_H

#include <vector>

#include "ElementErrorCode.h"

namespace MeshLib {
	class Mesh;

/**
 * \brief A collection of methods for testing mesh quality and correctness
 */
class MeshQualityController
{
public:
	MeshQualityController(MeshLib::Mesh &mesh);
	~MeshQualityController() {}

	/// Find mesh nodes that are not part of any element
	static std::vector<std::size_t> findUnusedMeshNodes(const MeshLib::Mesh &mesh);

	/// Removes nodes from the mesh that are not part of any element.
	/// Note: This method calls MeshQualityConroller::findUnusedMeshNodes() internally.
	static void removeUnusedMeshNodes(MeshLib::Mesh &mesh);

	/// Tests if elements are geometrically correct
	static std::vector<ElementErrorCode> testElementGeometry(const MeshLib::Mesh &mesh);

	/// Detailed output which ElementID is associated with which error(s)
	static void ElementErrorCodeOutput(const std::vector<ElementErrorCode> &error_codes);

private:

};

} // end namespace MeshLib

#endif //MESHQUALITYCONTROLLER_H
