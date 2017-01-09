/**
 * @copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
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
    BoundaryElementsAtPoint(MeshLib::Mesh const& mesh,
                           MeshNodeSearcher& mshNodeSearcher,
                           GeoLib::Point const& point);

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
