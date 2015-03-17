/**
 * \file
 * \author Karsten Rink
 * \date   2010-05-17
 * \brief  Implementation of the MshItem class.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MshItem.h"

/**
 * Constructor.
 * \param data The data associated with each column
 * \param parent The parent item in the tree
 * \param grid The mesh associated with this item
 */
MshItem::MshItem(const QList<QVariant> &data, TreeItem* parent, const MeshLib::Mesh* mesh)
: TreeItem(data, parent)
{
	_mesh_source = InSituLib::VtkMappedMeshSource::New();
	static_cast<InSituLib::VtkMappedMeshSource*>(_mesh_source)->SetMesh(mesh);
}

MshItem::~MshItem()
{
	_mesh_source->Delete();
}
