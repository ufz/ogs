/**
 * \file
 * \author Karsten Rink
 * \date   2010-05-17
 * \brief  Implementation of the MeshItem class.
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshItem.h"

#include "MeshLib/Vtk/VtkMappedMeshSource.h"

MeshItem::MeshItem(const QList<QVariant>& data, TreeItem* parent,
                   const MeshLib::Mesh* mesh)
    : TreeItem(data, parent)
{
    _mesh_source = MeshLib::VtkMappedMeshSource::New();
    static_cast<MeshLib::VtkMappedMeshSource*>(_mesh_source)->SetMesh(mesh);
}

MeshItem::~MeshItem()
{
    _mesh_source->Delete();
}
