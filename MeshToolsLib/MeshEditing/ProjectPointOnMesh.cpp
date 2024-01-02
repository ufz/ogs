/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
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

namespace MeshToolsLib
{
namespace ProjectPointOnMesh
{
MeshLib::Element const* getProjectedElement(
    std::vector<const MeshLib::Element*> const& elements,
    MathLib::Point3d const& node)
{
    auto is_right_of = [&node](MeshLib::Node const& a, MeshLib::Node const& b)
    { return GeoLib::getOrientation(node, a, b) == GeoLib::Orientation::CW; };

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

double getElevation(MeshLib::Element const& element,
                    MathLib::Point3d const& node)
{
    // mathematical description of the plane spanned by the 2d element
    // compute coefficients of the plane equation (Hesse normal form)
    // d = scalar_product(normal, element_node[0])
    auto const n = MeshLib::FaceRule::getSurfaceNormal(element).normalized();
    auto const d = n.dot(element.getNode(0)->asEigenVector3d());
    // insert node[0] and node[1] into plane equation and transpose the equation
    // to node[2]
    return (d - (node[0] * n[0] + node[1] * n[1])) / n[2];
}

}  // namespace ProjectPointOnMesh

}  // namespace MeshToolsLib
