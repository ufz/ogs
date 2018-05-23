/**
 * @copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
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
                            MeshNodeSearcher const& mshNodeSearcher,
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

    /// \copybrief BoundaryElementsSearcher::getBulkIDs()
    std::vector<std::pair<std::size_t, unsigned>> const& getBulkIDs() const;

private:
    MeshLib::Mesh const& _mesh;
    GeoLib::Point const& _point;
    std::vector<MeshLib::Element*> _boundary_elements;
    /// a vector of id pairs. The first item of each pair is the bulk element id
    /// and the second item is the face id of the corresponding bulk element.
    std::vector<std::pair<std::size_t, unsigned>> _bulk_ids;
};

}  // end namespace MeshGeoToolsLib
