/**
 * \file
 * \author Karsten Rink
 * \date   2010-05-17
 * \brief  Definition of the MshItem class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "TreeItem.h"

#include "MeshLib/Vtk/VtkMappedMeshSource.h"

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
    ///
    /// \param data The data associated with each column
    /// \param parent The parent item in the tree
    /// \param mesh The mesh associated with this item
    MshItem(const QList<QVariant>& data, TreeItem* parent,
            const MeshLib::Mesh* mesh);
    ~MshItem();

    /// Returns the mesh.
    MeshLib::Mesh const* getMesh() const { return _mesh_source->GetMesh(); }
    /// Returns the VTK object.
    MeshLib::VtkMappedMeshSource*  vtkSource() const { return _mesh_source; }

private:
    MeshLib::VtkMappedMeshSource * _mesh_source;
};
