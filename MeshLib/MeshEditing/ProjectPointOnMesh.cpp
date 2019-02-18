/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
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

/// Returns the element in which the given node is located when
/// projected onto a mesh, or nullptr if no such element was found.
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

/// Returns the z-coordinate of a point projected onto the plane defined
/// by a mesh element.
double getElevation(MeshLib::Element const& element,
                    MeshLib::Node const& node)
{
    MathLib::Vector3 const v = node - *element.getNode(0);
    MathLib::Vector3 const n =
        MeshLib::FaceRule::getSurfaceNormal(&element).getNormalizedVector();
    return node[2] - scalarProduct(n, v) * n[2];
}

void project(MeshLib::Mesh const& mesh,
             std::vector<MeshLib::Node*> const& nodes,
             double const default_value = 0)
{
    MeshLib::MeshElementGrid const grid(mesh);
    double const max_edge(mesh.getMaxEdgeLength());

    for (MeshLib::Node* node : nodes)
    {
        MathLib::Point3d min_vol{{(*node)[0] - max_edge, (*node)[1] - max_edge,
                                  -std::numeric_limits<double>::max()}};
        MathLib::Point3d max_vol{{(*node)[0] + max_edge, (*node)[1] + max_edge,
                                  std::numeric_limits<double>::max()}};
        std::vector<const MeshLib::Element*> const& elems =
            grid.getElementsInVolume(min_vol, max_vol);
        auto const* element = getProjectedElement(elems, *node);
        // centre of the pixel is located within a mesh element
        (*node)[2] =
            (element != nullptr) ? getElevation(*element, *node) : default_value;
    }
}

}  // namespace ProjectPointOnMesh

}  // end namespace MeshLib
