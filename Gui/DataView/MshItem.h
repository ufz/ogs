/**
 * \file MshItem.h
 * 17/05/2010 KR Initial implementation
 */

#ifndef MSHITEM_H
#define MSHITEM_H

#include "TreeItem.h"
#include "VtkMeshSource.h"

class VtkMeshSource;

namespace MeshLib {
	class Mesh;
}

/**
 * \brief A TreeItem containing a mesh and the associated vtk object used in the Mesh Model.
 * \sa TreeItem
 */
class MshItem : public TreeItem
{
public:
	/// Constructor, automatically generates VTK object of the given mesh.
	MshItem(const QList<QVariant> &data, TreeItem* parent, const MeshLib::Mesh* grid);
	~MshItem();

	/// Returns the mesh as a GridAdapter.
	const MeshLib::Mesh* getMesh() const { return this->_meshSource->GetMesh(); }
	/// Returns the VTK object.
	VtkMeshSource* vtkSource() const { return _meshSource; }

private:
	VtkMeshSource* _meshSource;
};

#endif //MSHITEM_H
