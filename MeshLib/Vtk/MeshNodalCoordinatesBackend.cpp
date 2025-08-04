/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

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
