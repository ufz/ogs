/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "GeoLib/AnalyticalGeometry.h"
#include "MathLib/Point3d.h"

#include "MeshLib/Elements/Quad.h"
#include "MeshLib/Elements/Tri.h"
#include "MeshLib/MeshSearch/MeshElementGrid.h"

namespace MeshLib
{
namespace ProjectPointOnMesh
{

MeshLib::Element const* getProjectedElement(
    std::vector<const MeshLib::Element*> const& elements,
    MeshLib::Node const& node)
{
    auto is_right_of = [&node](MeshLib::Node const& a, MeshLib::Node const& b) {
        return GeoLib::getOrientationFast(node, a, b) ==
               GeoLib::Orientation::CW;
    };

    for (auto const* e : elements)
    {
        auto const* nodes = e->getNodes();
        if (e->getGeomType() == MeshLib::MeshElemType::TRIANGLE)
        {
            auto const& a = *nodes[0];
            auto const& b = *nodes[1];
            auto const& c = *nodes[2];
            if (!is_right_of(a, b) && !is_right_of(b, c) && !is_right_of(c, a))
            {
                return e;
            }
        }
        else if (e->getGeomType() == MeshLib::MeshElemType::QUAD)
        {
            auto const& a = *nodes[0];
            auto const& b = *nodes[1];
            auto const& c = *nodes[2];
            auto const& d = *nodes[3];
            if (!is_right_of(a, b) && !is_right_of(b, c) &&
                !is_right_of(c, d) && !is_right_of(d, a))
            {
                return e;
            }
        }
    }
    return nullptr;
}

double getElevation(MeshLib::Element const& element, MeshLib::Node const& node)
{
    Eigen::Vector3d const v =
        Eigen::Map<Eigen::Vector3d const>(node.getCoords()) -
        Eigen::Map<Eigen::Vector3d const>(element.getNode(0)->getCoords());
    auto const n = MeshLib::FaceRule::getSurfaceNormal(&element).normalized();
    return node[2] - n.dot(v) * n[2];
}

}  // namespace ProjectPointOnMesh

}  // end namespace MeshLib
