// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "AngleSkewMetric.h"

#include <cmath>
#include <numbers>

#include "MathLib/MathTools.h"
#include "MeshLib/Node.h"

namespace MeshToolsLib
{
namespace
{
template <unsigned long N>
std::tuple<double, double> getMinMaxAngle(
    std::array<MeshLib::Node, N> const& nodes)
{
    double min_angle(2 * std::numbers::pi);
    double max_angle(0.0);

    for (decltype(N) i = 0; i < N; ++i)
    {
        const double angle(MathLib::getAngle(nodes[i], nodes[(i + 1) % N],
                                             nodes[(i + 2) % N]));
        min_angle = std::min(angle, min_angle);
        max_angle = std::max(angle, max_angle);
    }
    return {min_angle, max_angle};
}

double checkTriangle(MeshLib::Element const& elem)
{
    using namespace std::numbers;
    std::array const nodes = {*elem.getNode(0), *elem.getNode(1),
                              *elem.getNode(2)};
    auto const& [min_angle, max_angle] = getMinMaxAngle(nodes);
    return std::max((max_angle - pi / 3) / 2, (pi / 3 - min_angle)) * 3 / pi;
}

double checkQuad(MeshLib::Element const& elem)
{
    std::array const nodes = {*elem.getNode(0), *elem.getNode(1),
                              *elem.getNode(2), *elem.getNode(3)};
    auto const& [min_angle, max_angle] = getMinMaxAngle(nodes);

    using namespace std::numbers;
    return std::max((max_angle - pi / 2), (pi / 2 - min_angle)) * 2 / pi;
}

double checkTetrahedron(MeshLib::Element const& elem)
{
    std::array<double, 4> min;
    std::array<double, 4> max;
    for (auto face_number = 0; face_number < 4; ++face_number)
    {
        std::unique_ptr<MeshLib::Element const> face{elem.getFace(face_number)};
        std::array const nodes = {*face->getNode(0), *face->getNode(1),
                                  *face->getNode(2)};
        std::tie(min[face_number], max[face_number]) = getMinMaxAngle(nodes);
    }

    double const min_angle = *std::min_element(min.begin(), min.end());
    double const max_angle = *std::max_element(max.begin(), max.end());

    using namespace std::numbers;
    return std::max((max_angle - pi / 3) / 2, (pi / 3 - min_angle)) * 3 / pi;
}

double checkHexahedron(MeshLib::Element const& elem)
{
    std::array<double, 6> min;
    std::array<double, 6> max;
    for (auto face_number = 0; face_number < 6; ++face_number)
    {
        std::unique_ptr<MeshLib::Element const> face{elem.getFace(face_number)};
        std::array const nodes = {*face->getNode(0), *face->getNode(1),
                                  *face->getNode(2), *face->getNode(3)};
        std::tie(min[face_number], max[face_number]) = getMinMaxAngle(nodes);
    }

    double const min_angle = *std::min_element(min.begin(), min.end());
    double const max_angle = *std::max_element(max.begin(), max.end());

    using namespace std::numbers;
    return std::max((max_angle - pi / 2), (pi / 2 - min_angle)) * 2 / pi;
}

double checkPrism(MeshLib::Element const& elem)
{
    // face 0: triangle (0,1,2)
    std::unique_ptr<MeshLib::Element const> f0{elem.getFace(0)};
    std::array const nodes_f0 = {*f0->getNode(0), *f0->getNode(1),
                                 *f0->getNode(2)};
    auto const& [min_angle_tri0, max_angle_tri0] = getMinMaxAngle(nodes_f0);

    // face 4: triangle (3,4,5)
    std::unique_ptr<MeshLib::Element const> f4{elem.getFace(4)};
    std::array const nodes_f4 = {*f4->getNode(0), *f4->getNode(1),
                                 *f4->getNode(2)};
    auto const& [min_angle_tri1, max_angle_tri1] = getMinMaxAngle(nodes_f4);

    auto const min_angle_tri = std::min(min_angle_tri0, min_angle_tri1);
    auto const max_angle_tri = std::max(max_angle_tri0, max_angle_tri1);

    using namespace std::numbers;
    double const tri_criterion =
        std::max((max_angle_tri - pi / 3) / 2, (pi / 3 - min_angle_tri)) * 3 /
        pi;

    std::array<double, 3> min;
    std::array<double, 3> max;
    for (int i = 1; i < 4; ++i)
    {
        std::unique_ptr<MeshLib::Element const> f{elem.getFace(i)};
        std::array const nodes = {*f->getNode(0), *f->getNode(1),
                                  *f->getNode(2), *f->getNode(3)};
        std::tie(min[i - 1], max[i - 1]) = getMinMaxAngle(nodes);
    }

    double const min_angle_quad = *std::min_element(min.begin(), min.end());
    double const max_angle_quad = *std::max_element(max.begin(), max.end());

    using namespace std::numbers;
    double const quad_criterion =
        std::max((max_angle_quad - pi / 2), (pi / 2 - min_angle_quad)) * 2 / pi;

    return std::min(tri_criterion, quad_criterion);
}

}  // end unnamed namespace

void AngleSkewMetric::calculateQuality()
{
    for (auto const e : _mesh.getElements())
    {
        switch (e->getGeomType())
        {
            case MeshLib::MeshElemType::LINE:
                _element_quality_metric[e->getID()] = -1.0;
                break;
            case MeshLib::MeshElemType::TRIANGLE:
                _element_quality_metric[e->getID()] = checkTriangle(*e);
                break;
            case MeshLib::MeshElemType::QUAD:
                _element_quality_metric[e->getID()] = checkQuad(*e);
                break;
            case MeshLib::MeshElemType::TETRAHEDRON:
                _element_quality_metric[e->getID()] = checkTetrahedron(*e);
                break;
            case MeshLib::MeshElemType::HEXAHEDRON:
                _element_quality_metric[e->getID()] = checkHexahedron(*e);
                break;
            case MeshLib::MeshElemType::PRISM:
                _element_quality_metric[e->getID()] = checkPrism(*e);
                break;
            default:
                break;
        }
    }
}

}  // namespace MeshToolsLib
