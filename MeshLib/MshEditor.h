/**
 * \file MshEditor.h
 * 2011/06/15 KR Initial implementation
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
