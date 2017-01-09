/**
 * \file
 * \author Karsten Rink
 * \date   2010-05-17
 * \brief  Implementation of the MshItem class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MshItem.h"
#include "MeshLib/Vtk/VtkMappedMeshSource.h"

MshItem::MshItem(const QList<QVariant>& data, TreeItem* parent,
                 const MeshLib::Mesh* mesh)
    : TreeItem(data, parent)
{
    _mesh_source = MeshLib::VtkMappedMeshSource::New();
    static_cast<MeshLib::VtkMappedMeshSource*>(_mesh_source)->SetMesh(mesh);
}

MshItem::~MshItem()
{
    _mesh_source->Delete();
}
