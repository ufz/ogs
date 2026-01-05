// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
