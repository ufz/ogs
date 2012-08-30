/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 * \file MshItem.cpp
 *
 * Created on 2010-05-17 by Karsten Rink
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
