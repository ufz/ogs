/**
 * \file
 * \author Lars Bilke
 * \date   2009-10-19
 * \brief  Definition of the MshModel class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>

#include "MeshLib/MeshEnums.h"

#include "TreeModel.h"

namespace DataHolderLib
{
class Project;
}

namespace MeshLib {
class Mesh;
}

class vtkUnstructuredGridAlgorithm;

/**
 * The MshModel is a Qt model which represents Mesh objects.
 */
class MshModel : public TreeModel
{
    Q_OBJECT

public:
    MshModel(DataHolderLib::Project& project, QObject* parent = nullptr);

    /// Adds a new mesh
    void addMesh(std::unique_ptr<MeshLib::Mesh> mesh);

    /// Returns the number of columns used for the data list
    int columnCount(const QModelIndex& parent = QModelIndex()) const override;

public slots:
    /// Adds a new mesh (using no unique_ptr as this interferes with Qt's signal/slot policy)
    void addMesh(MeshLib::Mesh* mesh);
    /// Returns the mesh with the given index.
    const MeshLib::Mesh* getMesh(const QModelIndex &idx) const;
    /// Returns the mesh with the given name.
    const MeshLib::Mesh* getMesh(const std::string &name) const;
    /// Removes the mesh with the given index.
    bool removeMesh(const QModelIndex &idx);
    /// Removes the mesh with the given name.
    bool removeMesh(const std::string &name);
    /// Updates the model/view for a mesh.
    void updateMesh(MeshLib::Mesh*);
    /// Updates the model based on the ProjectData-object
    void updateModel();
    /// Returns the VTK source item for the mesh with the given index.
    vtkUnstructuredGridAlgorithm* vtkSource(const QModelIndex &idx) const;
    /// Returns the VTK source item for the mesh with the given name.
    vtkUnstructuredGridAlgorithm* vtkSource(const std::string &name) const;

private:
    /// Adds the mesh to the GUI-Mesh-Model und -View
    void addMeshObject(const MeshLib::Mesh* mesh);

    /// Checks if the name of the mesh is already exists, if so it generates a unique name.
    //bool isUniqueMeshName(std::string &name);
    DataHolderLib::Project& _project;

    /// Creates a static map of all element type name-strings in QVariant format
    static std::map<MeshLib::MeshElemType, QVariant> createMeshElemTypeMap();

    static const std::map<MeshLib::MeshElemType, QVariant> elem_type_map;
    static const QVariant element_str;

signals:
    void meshAdded(MshModel*, const QModelIndex&);
    void meshRemoved(MshModel*, const QModelIndex&);
};
