// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "MeshNodalCoordinatesBackend.h"

#include <algorithm>

#include "MeshLib/Node.h"

namespace MeshLib
{
MeshNodalCoordinatesBackend::MeshNodalCoordinatesBackend(
    std::vector<Node*> const& nodes)
    : nodes_(nodes)
{
}

double MeshNodalCoordinatesBackend::map(vtkIdType idx) const
{
    return nodes_[idx / 3]->operator[](idx % 3);
}

void MeshNodalCoordinatesBackend::mapTuple(vtkIdType tupleId,
                                           double* tuple) const
{
    std::copy_n(nodes_[tupleId]->data(), 3, tuple);
}
}  // namespace MeshLib
