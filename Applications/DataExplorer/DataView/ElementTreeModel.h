// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <array>
#include "Base/TreeModel.h"

class vtkUnstructuredGridAlgorithm;

namespace MeshLib {
    class Mesh;
}

/**
 * \brief A model for the display of information concerning element information implemented as a TreeModel.
 * \sa TreeModel, ElementTreeView, TreeItem
 */
class ElementTreeModel : public TreeModel
{
    Q_OBJECT

public:
    explicit ElementTreeModel(QObject* parent = nullptr);
    ~ElementTreeModel() override;

    vtkUnstructuredGridAlgorithm const* getSource() const { return _mesh_source; };

public slots:
    /// Clears the tree.
    void clearView();

    /// Displays information of the element with the given index from the given grid.
    void setElement(vtkUnstructuredGridAlgorithm const*const grid, const unsigned elem_index);

    /// Displays information of the given mesh.
    void setMesh(MeshLib::Mesh const& mesh);

private:
    vtkUnstructuredGridAlgorithm const* _mesh_source{nullptr};
};
