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
    // *** rotate polygon to xy_plane
    Eigen::Vector3d normal;
    GeoLib::Polygon rot_polygon(GeoLib::rotatePolygonToXY(polygon, normal));

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

    // *** mark rotated nodes
    std::vector<bool> outside(rotated_nodes.size(), true);
    for (std::size_t k(0); k < rotated_nodes.size(); k++)
    {
        if (rot_polygon.isPntInPolygon(*(rotated_nodes[k])))
        {
            outside[k] = false;
        }
    }

    for (auto& rotated_node : rotated_nodes)
    {
        delete rotated_node;
    }

    std::vector<GeoLib::Point*>& rot_polygon_pnts(
        const_cast<std::vector<GeoLib::Point*>&>(rot_polygon.getPointsVec()));
    for (auto& rot_polygon_pnt : rot_polygon_pnts)
    {
        delete rot_polygon_pnt;
    }

    return outside;
}

}  // end namespace MeshGeoToolsLib
