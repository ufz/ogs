/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 * \file MshEditor.h
 *
 * Created on 2011-06-15 by Karsten Rink
 */

#ifndef MSHEDITOR_H
#define MSHEDITOR_H

#include <cstddef>
#include <vector>

namespace GeoLib {
	class PointWithID;
}

namespace MeshLib
{
class Mesh;
}

class GridAdapter;

/**
 * \brief A set of tools for manipulating existing meshes
 */
class MshEditor
{
public:
	MshEditor() {}
	~MshEditor() {}

	static MeshLib::Mesh* removeMeshNodes(MeshLib::Mesh* mesh,
	                                         const std::vector<size_t> &nodes);

	static std::vector<GeoLib::PointWithID*> getSurfaceNodes(const MeshLib::Mesh &mesh);

	static MeshLib::Mesh* getMeshSurface(const MeshLib::Mesh &mesh);

private:
};

#endif //MSHEDITOR_H
