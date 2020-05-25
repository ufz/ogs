/**
 * \file
 * \author Karsten Rink
 * \date   2010-05-17
 * \brief  Implementation of the MeshItem class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
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
    mesh_source_ = MeshLib::VtkMappedMeshSource::New();
    static_cast<MeshLib::VtkMappedMeshSource*>(mesh_source_)->SetMesh(mesh);
}

MeshItem::~MeshItem()
{
    mesh_source_->Delete();
}
