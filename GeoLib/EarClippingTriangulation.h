/**
 * \file
 * \author Thomas Fischer
 * \date   2011-02-23
 * \brief  Definition of the EarClippingTriangulation class.
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <list>
#include <vector>

#include "AnalyticalGeometry.h"

namespace GeoLib
{
class Polygon;
class Triangle;

class EarClippingTriangulation final
{
public:
    EarClippingTriangulation(GeoLib::Polygon const& polygon,
                             std::list<GeoLib::Triangle>& triangles,
                             bool rot = true);
    ~EarClippingTriangulation();

private:
    /**
     * copies the points of the polygon to the vector _pnts
     */
    inline void copyPolygonPoints(GeoLib::Polygon const& polygon);
    inline void ensureCWOrientation();

    inline bool isEar(std::size_t v0, std::size_t v1, std::size_t v2) const;

    inline void initVertexList();
    inline void initLists();
    inline void clipEars();
    void addLastTriangle();

    /**
     * a copy of the polygon points
     */
    std::vector<GeoLib::Point*> _pnts;
    std::list<std::size_t> _vertex_list;
    std::list<std::size_t> _convex_vertex_list;
    std::list<std::size_t> _ear_list;

    /**
     * triangles of the triangulation (maybe in the wrong orientation)
     */
    std::list<GeoLib::Triangle> _triangles;

    GeoLib::Orientation _original_orientation;
};
}  // end namespace GeoLib
