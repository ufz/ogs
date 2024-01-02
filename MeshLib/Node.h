/**
 * \file
 * \author Karsten Rink
 * \date   2012-05-02
 * \brief  Definition of the Node class.
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <cstdlib>
#include <limits>
#include <vector>

#include "MathLib/Point3dWithID.h"

namespace ApplicationUtils
{
    class NodeWiseMeshPartitioner;
}

namespace MeshToolsLib
{
class MeshRevision;
}

namespace MeshLib
{
/**
 * A mesh node with coordinates in 3D space.
 */
class Node : public MathLib::Point3dWithID
{
    /* friend classes: */
    friend class Mesh;
    friend class MeshToolsLib::MeshRevision;

public:
    /// Constructor using a coordinate array
    explicit Node(const double coords[3],
                  std::size_t id = std::numeric_limits<std::size_t>::max());

    /// Constructor using a coordinate array
    explicit Node(std::array<double, 3> const& coords,
                  std::size_t id = std::numeric_limits<std::size_t>::max());

    /// Constructor using single coordinates
    Node(double x, double y, double z, std::size_t id = std::numeric_limits<std::size_t>::max());

    /// Copy constructor
    Node(const Node &node);

}; /* class */
}  // namespace MeshLib
