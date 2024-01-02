/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */
#pragma once

#include <vector>

namespace GeoLib
{
class Polyline;
}

namespace MeshLib
{
class Mesh;
class Element;
}  // namespace MeshLib

namespace MeshGeoToolsLib
{
class MeshNodeSearcher;

/**
 * This class collects element edges located along a polyline.
 * Note that internal edges are not collected in this class.
 */
class BoundaryElementsAlongPolyline
{
public:
    /**
     * Constructor
     * @param mesh             a mesh object
     * @param mshNodeSearcher  a MeshNodeSearcher object which is internally
     * used to search mesh nodes
     * @param ply              a polyline object where edges are searched
     */
    BoundaryElementsAlongPolyline(MeshLib::Mesh const& mesh,
                                  MeshNodeSearcher const& mshNodeSearcher,
                                  GeoLib::Polyline const& ply);

    /// destructor
    virtual ~BoundaryElementsAlongPolyline();

    /**
     * Deploying this method the user can get access to the underlying
     * GeoLib::Polyline.
     * @return the underlying GeoLib::Polyline
     */
    GeoLib::Polyline const& getPolyline() const { return _ply; }

    /**
     * Return the vector of boundary elements (i.e. edges). The elements are
     * sorted according to their distance to the starting point of the given
     * polyline.
     */
    std::vector<MeshLib::Element*> const& getBoundaryElements() const
    {
        return _boundary_elements;
    }

private:
    GeoLib::Polyline const& _ply;
    std::vector<MeshLib::Element*> _boundary_elements;
};

}  // end namespace MeshGeoToolsLib
