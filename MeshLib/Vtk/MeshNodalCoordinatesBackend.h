/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */
#pragma once

#include <vtkType.h>

#include <vector>

namespace MeshLib
{
class Node;
}

namespace MeshLib
{
struct MeshNodalCoordinatesBackend
{
    explicit MeshNodalCoordinatesBackend(std::vector<Node*> const& nodes);

    double map(vtkIdType idx) const;
    void mapTuple(vtkIdType tupleId, double* tuple) const;

private:
    std::vector<Node*> const& nodes_;
};
}  // namespace MeshLib
