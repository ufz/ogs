/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <vector>

namespace GeoLib
{
class Point;
}

namespace MeshLib
{
class Mesh;
class Element;
}

namespace MeshGeoToolsLib
{
class MeshNodeSearcher;

/// This class collects point elements located at a given point elements.
class BoundaryElementsAtPoint final
{
public:
    /// Constructor
    /// \param mesh             a mesh object
    /// \param mshNodeSearcher  a MeshNodeSearcher object which is internally
    /// used to search mesh nodes
    /// \param point            a point object where edges are searched
    /// \param multiple_nodes_allowed Allows to find multiple nodes within the
    /// search radius, the nearest node is returned (as point element). This
    /// enables to specify larger search radius to find possible other
    /// geometries that don't match exactly to the mesh.
    BoundaryElementsAtPoint(MeshLib::Mesh const& mesh,
                            MeshNodeSearcher const& mshNodeSearcher,
                            GeoLib::Point const& point,
                            const bool multiple_nodes_allowed);

    ~BoundaryElementsAtPoint();

    MeshLib::Mesh const& getMesh() const
    {
        return _mesh;
    }

    GeoLib::Point const& getPoint() const
    {
        return _point;
    }

    /// Return the vector of boundary elements (i.e. points).
    std::vector<MeshLib::Element*> const& getBoundaryElements() const
    {
        return _boundary_elements;
    }

private:
    MeshLib::Mesh const& _mesh;
    GeoLib::Point const& _point;
    std::vector<MeshLib::Element*> _boundary_elements;
};

}  // end namespace MeshGeoToolsLib
