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

	/// Removes nodes from the mesh that are not part of any element.
	static void removeUnusedMeshNodes(MeshLib::Mesh &mesh);

	/// Tests if elements are geometrically correct
	static void testElementGeometry(MeshLib::Mesh &mesh);

private:

};

} // end namespace MeshLib

#endif //MESHQUALITYCONTROLLER_H
