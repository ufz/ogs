/**
 * \file   AngleSkewMetric.cpp
 * \author Thomas Fischer
 * \date   2011-03-17
 * \brief  Implementation of the AngleSkewMetric class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
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
            _element_quality_metric[k] = checkTriangle (elem);
            break;
        case MeshElemType::QUAD:
            _element_quality_metric[k] = checkQuad (elem);
            break;
        case MeshElemType::TETRAHEDRON:
            _element_quality_metric[k] = checkTetrahedron (elem);
            break;
        case MeshElemType::HEXAHEDRON:
            _element_quality_metric[k] = checkHexahedron (elem);
            break;
        case MeshElemType::PRISM:
            _element_quality_metric[k] = checkPrism (elem);
            break;
        default:
            break;
        }
    }
}

double AngleSkewMetric::checkTriangle (Element const& elem) const
{
    double const* const node0 (elem.getNode(0)->getCoords());
    double const* const node1 (elem.getNode(1)->getCoords());
    double const* const node2 (elem.getNode(2)->getCoords());

    double min_angle (two_pi), max_angle (0.0);
    getMinMaxAngleFromTriangle (node0, node1, node2, min_angle, max_angle);

    return 1.0 -
           std::max((max_angle - third_pi) / two_thirds_pi,
                    (third_pi - min_angle) / third_pi);
}

double AngleSkewMetric::checkQuad (Element const& elem) const
{
    double const* const node0 (elem.getNode(0)->getCoords());
    double const* const node1 (elem.getNode(1)->getCoords());
    double const* const node2 (elem.getNode(2)->getCoords());
    double const* const node3 (elem.getNode(3)->getCoords());

    double min_angle (two_pi);
    double max_angle (0.0);

    getMinMaxAngleFromQuad (node0, node1, node2, node3, min_angle, max_angle);

    return 1.0 -
           std::max((max_angle - two_pi) / (-pi), (two_pi - min_angle) / (two_pi));
}

double AngleSkewMetric::checkTetrahedron (Element const& elem) const
{
    double const* const node0 (elem.getNode(0)->getCoords());
    double const* const node1 (elem.getNode(1)->getCoords());
    double const* const node2 (elem.getNode(2)->getCoords());
    double const* const node3 (elem.getNode(3)->getCoords());

    double min_angle (two_pi);
    double max_angle (0.0);

    // first triangle (0,1,2)
    getMinMaxAngleFromTriangle(node0, node1, node2, min_angle, max_angle);
    // second triangle (0,1,3)
    getMinMaxAngleFromTriangle(node0, node1, node3, min_angle, max_angle);
    // third triangle (0,2,3)
    getMinMaxAngleFromTriangle(node0, node2, node3, min_angle, max_angle);
    // fourth triangle (1,2,3)
    getMinMaxAngleFromTriangle(node1, node2, node3, min_angle, max_angle);

    return 1.0 - std::max((max_angle - two_pi) / two_thirds_pi,
                          (third_pi - min_angle) / third_pi);
}

double AngleSkewMetric::checkHexahedron (Element const& elem) const
{
    double const* const node0 (elem.getNode(0)->getCoords());
    double const* const node1 (elem.getNode(1)->getCoords());
    double const* const node2 (elem.getNode(2)->getCoords());
    double const* const node3 (elem.getNode(3)->getCoords());
    double const* const node4 (elem.getNode(4)->getCoords());
    double const* const node5 (elem.getNode(5)->getCoords());
    double const* const node6 (elem.getNode(6)->getCoords());
    double const* const node7 (elem.getNode(7)->getCoords());

    double min_angle (two_pi);
    double max_angle (0.0);

    // first surface (0,1,2,3)
    getMinMaxAngleFromQuad (node0, node1, node2, node3, min_angle, max_angle);
    // second surface (0,3,7,4)
    getMinMaxAngleFromQuad (node0, node3, node7, node4, min_angle, max_angle);
    // third surface (4,5,6,7)
    getMinMaxAngleFromQuad (node4, node5, node6, node7, min_angle, max_angle);
    // fourth surface (5,1,2,6)
    getMinMaxAngleFromQuad (node5, node1, node2, node6, min_angle, max_angle);
    // fifth surface (5,1,0,4)
    getMinMaxAngleFromQuad (node5, node1, node0, node4, min_angle, max_angle);
    // sixth surface (6,2,3,7)
    getMinMaxAngleFromQuad (node6, node2, node3, node7, min_angle, max_angle);

    return 1.0 -
           std::max((max_angle - two_pi) / (-pi), (two_pi - min_angle) / two_pi);
}

double AngleSkewMetric::checkPrism (Element const& elem) const
{
    double const* const node0 (elem.getNode(0)->getCoords());
    double const* const node1 (elem.getNode(1)->getCoords());
    double const* const node2 (elem.getNode(2)->getCoords());
    double const* const node3 (elem.getNode(3)->getCoords());
    double const* const node4 (elem.getNode(4)->getCoords());
    double const* const node5 (elem.getNode(5)->getCoords());

    double min_angle_tri (two_pi);
    double max_angle_tri (0.0);

    // first triangle (0,1,2)
    getMinMaxAngleFromTriangle (node0, node1, node2, min_angle_tri, max_angle_tri);
    // second surface (3,4,5)
    getMinMaxAngleFromTriangle (node3, node4, node5, min_angle_tri, max_angle_tri);

    double tri_criterion (1.0 - std::max((max_angle_tri - two_pi) / two_thirds_pi,
                                         (third_pi - min_angle_tri) / third_pi));

    double min_angle_quad (two_pi);
    double max_angle_quad (0.0);
    // surface (0,3,4,1)
    getMinMaxAngleFromQuad (node0, node3, node4, node1, min_angle_quad, max_angle_quad);
    // surface (2,5,3,0)
    getMinMaxAngleFromQuad (node2, node5, node3, node0, min_angle_quad, max_angle_quad);
    // surface (1,2,5,4)
    getMinMaxAngleFromQuad (node1, node2, node5, node4, min_angle_quad, max_angle_quad);

    double quad_criterion (1.0 - std::max((max_angle_quad - two_pi) / (-pi),
                                          (two_pi - min_angle_quad) / two_pi));

    return std::min (tri_criterion, quad_criterion);
}

void AngleSkewMetric::getMinMaxAngleFromQuad (
        double const* const n0, double const* const n1,
        double const* const n2, double const* const n3,
        double &min_angle, double &max_angle) const
{
    const double* nodes[4] = {n0, n1, n2, n3};
    for (unsigned i=0; i<4; ++i)
    {
        const double angle (MathLib::getAngle (nodes[i], nodes[(i+1)%4], nodes[(i+2)%4]));
        min_angle = std::min(angle, min_angle);
        max_angle = std::max(angle, max_angle);
    }
}

void AngleSkewMetric::getMinMaxAngleFromTriangle(double const* const n0,
                                                          double const* const n1,
                                                          double const* const n2,
                                                          double &min_angle,
                                                          double &max_angle) const
{
    const double* nodes[3] = {n0, n1, n2};
    for (unsigned i=0; i<3; ++i)
    {
        const double angle (MathLib::getAngle (nodes[i], nodes[(i+1)%3], nodes[(i+2)%3]));
        min_angle = std::min(angle, min_angle);
        max_angle = std::max(angle, max_angle);
    }
}

} // end namespace MeshLib
