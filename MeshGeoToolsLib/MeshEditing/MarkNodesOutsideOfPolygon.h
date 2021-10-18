/*
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <algorithm>
#include <vector>

#include "GeoLib/AnalyticalGeometry.h"
#include "GeoLib/Polygon.h"

#include "MeshLib/Node.h"

namespace MeshGeoToolsLib
{
std::vector<bool> markNodesOutSideOfPolygon(
    std::vector<MeshLib::Node*> const& nodes, GeoLib::Polygon const& polygon)
{
    // *** rotate polygon points to xy-plane
    auto [rotated_polygon_points, normal] =
        GeoLib::rotatePolygonPointsToXY(polygon);

    // *** rotate mesh nodes to xy-plane
    // 1 copy all mesh nodes to GeoLib::Points
    std::vector<GeoLib::Point*> rotated_nodes;
    for (auto node : nodes)
    {
        rotated_nodes.push_back(new GeoLib::Point(*node, node->getID()));
    }
    // 2 rotate the Points
    Eigen::Matrix3d const rot_mat = GeoLib::computeRotationMatrixToXY(normal);
    GeoLib::rotatePoints(rot_mat, rotated_nodes);
    // 3 set z coord to zero
    std::for_each(rotated_nodes.begin(), rotated_nodes.end(),
                  [](GeoLib::Point* p) { (*p)[2] = 0.0; });

    std::vector<bool> outside(rotated_nodes.size(), true);
    // *** mark rotated nodes inside rotated polygon
    {
        // create new polygon using the rotated points
        GeoLib::Polyline rotated_polyline(rotated_polygon_points);
        for (std::size_t k(0); k < polygon.getNumberOfPoints(); k++)
        {
            rotated_polyline.addPoint(k);
        }
        rotated_polyline.addPoint(0);
        GeoLib::Polygon const rotated_polygon(rotated_polyline);

        for (std::size_t k(0); k < rotated_nodes.size(); k++)
        {
            if (rotated_polygon.isPntInPolygon(*(rotated_nodes[k])))
            {
                outside[k] = false;
            }
        }
    }

    for (auto& rotated_node : rotated_nodes)
    {
        delete rotated_node;
    }

    for (auto& rot_polygon_pnt : rotated_polygon_points)
    {
        delete rot_polygon_pnt;
    }

    return outside;
}

}  // end namespace MeshGeoToolsLib
