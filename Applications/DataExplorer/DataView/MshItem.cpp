/**
 * \file
 * \author Karsten Rink
 * \date   2010-05-17
 * \brief  Implementation of the MshItem class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MshItem.h"
#include "VtkMeshSource.h"

/**
 * Constructor.
 * \param data The data associated with each column
 * \param parent The parent item in the tree
 * \param grid The mesh associated with this item
 */
MshItem::MshItem(const QList<QVariant> &data, TreeItem* parent, const MeshLib::Mesh* mesh)
	: TreeItem(data, parent)
{
	_meshSource = VtkMeshSource::New();
	_meshSource->SetMesh(mesh);
}

MshItem::~MshItem()
{
	_meshSource->Delete();
}
