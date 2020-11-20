/**
 * \file
 * \author Thomas Fischer
 * \date   2011-03-17
 * \brief  Implementation of the AngleSkewMetric class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "AngleSkewMetric.h"

#include <cmath>
#include <boost/math/constants/constants.hpp>

#include "MeshLib/Node.h"

#include "MathLib/MathTools.h"

using namespace boost::math::double_constants;

namespace MeshLib
{
AngleSkewMetric::AngleSkewMetric(Mesh const& mesh) :
    ElementQualityMetric(mesh)
{}

void AngleSkewMetric::calculateQuality ()
{
    const std::vector<MeshLib::Element*>& elements(_mesh.getElements());
    const std::size_t nElements (_mesh.getNumberOfElements());

    for (std::size_t k(0); k < nElements; k++)
    {
        Element const& elem (*elements[k]);
        switch (elem.getGeomType())
        {
        case MeshElemType::LINE:
            _element_quality_metric[k] = -1.0;
            break;
        case MeshElemType::TRIANGLE:
            _element_quality_metric[k] = checkTriangle(elem);
            break;
        case MeshElemType::QUAD:
            _element_quality_metric[k] = checkQuad(elem);
            break;
        case MeshElemType::TETRAHEDRON:
            _element_quality_metric[k] = checkTetrahedron(elem);
            break;
        case MeshElemType::HEXAHEDRON:
            _element_quality_metric[k] = checkHexahedron(elem);
            break;
        case MeshElemType::PRISM:
            _element_quality_metric[k] = checkPrism(elem);
            break;
        default:
            break;
        }
    }
}

double AngleSkewMetric::checkTriangle(Element const& elem) const
{
    auto const& [min_angle, max_angle] = getMinMaxAngleFromTriangle(
        *elem.getNode(0), *elem.getNode(1), *elem.getNode(2));
    return std::max((max_angle - third_pi) / two_thirds_pi,
                    (third_pi - min_angle) / third_pi);
}

double AngleSkewMetric::checkQuad(Element const& elem) const
{
    auto const& [min_angle, max_angle] = getMinMaxAngleFromQuad(
        *elem.getNode(0), *elem.getNode(1), *elem.getNode(2), *elem.getNode(3));


    return std::max((max_angle - half_pi) / half_pi,
                    (half_pi - min_angle) / half_pi);
}

double AngleSkewMetric::checkTetrahedron(Element const& elem) const
{
    std::array<double, 4> min;
    std::array<double, 4> max;
    for (auto face_number = 0; face_number < 4; ++face_number)
    {
        auto const& face = *elem.getFace(face_number);
        std::tie(min[face_number], max[face_number]) = getMinMaxAngleFromTri(
            *face.getNode(0), *face.getNode(1), *face.getNode(2));
    }

    double const min_angle = *std::min_element(min.begin(), min.end());
    double const max_angle = *std::max_element(max.begin(), max.end());

    return std::max((max_angle - third_pi) / two_thirds_pi,
                    (third_pi - min_angle) / third_pi);
}

double AngleSkewMetric::checkHexahedron(Element const& elem) const
{
    std::array<double, 6> min;
    std::array<double, 6> max;
    for (auto face_number = 0; face_number < 6; ++face_number)
    {
        auto const& face = *elem.getFace(face_number);
        std::tie(min[face_number], max[face_number]) =
            getMinMaxAngleFromQuad(*face.getNode(0), *face.getNode(1),
                                   *face.getNode(2), *face.getNode(3));
    }

    double const min_angle = *std::min_element(min.begin(), min.end());
    double const max_angle = *std::max_element(max.begin(), max.end());

    return std::max((max_angle - half_pi) / half_pi,
                    (half_pi - min_angle) / half_pi);
}

double AngleSkewMetric::checkPrism(Element const& elem) const
{
    // first triangle (0,1,2)
    auto const& [min_angle_tri0, max_angle_tri0] = getMinMaxAngleFromTriangle(
        *elem.getNode(0), *elem.getNode(1), *elem.getNode(2));
    // second surface (3,4,5)
    auto const& [min_angle_tri1, max_angle_tri1] = getMinMaxAngleFromTriangle(
        *elem.getNode(3), *elem.getNode(4), *elem.getNode(5));
    double const min_angle_tri = std::min(min_angle_tri0, min_angle_tri1);
    double const max_angle_tri = std::max(max_angle_tri0, max_angle_tri1);

    double const tri_criterion(
        std::max((max_angle_tri - third_pi) / two_thirds_pi,
                 (third_pi - min_angle_tri) / third_pi));

    // surface (0,3,4,1)
    auto const& [min_angle_quad0, max_angle_quad0] = getMinMaxAngleFromQuad(
        *elem.getNode(0), *elem.getNode(3), *elem.getNode(4), *elem.getNode(1));
    // surface (2,5,3,0)
    auto const& [min_angle_quad1, max_angle_quad1] = getMinMaxAngleFromQuad(
        *elem.getNode(2), *elem.getNode(5), *elem.getNode(3), *elem.getNode(0));
    // surface (1,2,5,4)
    auto const& [min_angle_quad2, max_angle_quad2] = getMinMaxAngleFromQuad(
        *elem.getNode(1), *elem.getNode(2), *elem.getNode(5), *elem.getNode(4));

    double const min_angle_quad =
        std::min({min_angle_quad0, min_angle_quad1, min_angle_quad2});
    double const max_angle_quad =
        std::max({max_angle_quad0, max_angle_quad1, max_angle_quad2});
}

std::tuple<double, double> AngleSkewMetric::getMinMaxAngleFromQuad(
    MeshLib::Node const& n0, MeshLib::Node const& n1, MeshLib::Node const& n2,
    MeshLib::Node const& n3) const
{
    double min_angle(two_pi);
    double max_angle(0.0);

    std::array const nodes = {n0, n1, n2, n3};
    for (unsigned i = 0; i < nodes.size(); ++i)
    {
        const double angle(MathLib::getAngle(nodes[i].getCoords(),
                                             nodes[(i + 1) % 4].getCoords(),
                                             nodes[(i + 2) % 4].getCoords()));
        min_angle = std::min(angle, min_angle);
        max_angle = std::max(angle, max_angle);
    }
    return {min_angle, max_angle};
}

std::tuple<double, double> AngleSkewMetric::getMinMaxAngleFromTriangle(
    MeshLib::Node const& n0, MeshLib::Node const& n1,
    MeshLib::Node const& n2) const
{
    double min_angle(two_pi);
    double max_angle(0.0);

    std::array nodes = {n0, n1, n2};
    for (unsigned i=0; i<3; ++i)
    {
        const double angle(MathLib::getAngle(nodes[i].getCoords(),
                                             nodes[(i + 1) % 3].getCoords(),
                                             nodes[(i + 2) % 3].getCoords()));
        min_angle = std::min(angle, min_angle);
        max_angle = std::max(angle, max_angle);
    }
    return {min_angle, max_angle};
}

} // end namespace MeshLib
