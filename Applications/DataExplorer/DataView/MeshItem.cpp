// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
