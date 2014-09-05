/**
 * \file
 * \author Karsten Rink
 * \date   2010-05-17
 * \brief  Definition of the MshItem class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
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

	/// Returns the mesh.
	const MeshLib::Mesh* getMesh() const { return this->_meshSource->GetMesh(); }
	/// Returns the VTK object.
	VtkMeshSource* vtkSource() const { return _meshSource; }

private:
	VtkMeshSource* _meshSource;
};

#endif //MSHITEM_H
