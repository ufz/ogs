/**
 * \file   EdgeRatioMetric.cpp
 * \author Thomas Fischer
 * \date   2011-03-03
 * \brief  Implementation of the EdgeRatioMetric class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "EdgeRatioMetric.h"

#include "MathLib/MathTools.h"
#include "MeshLib/Node.h"

namespace MeshLib
{
EdgeRatioMetric::EdgeRatioMetric(Mesh const& mesh) :
    ElementQualityMetric(mesh)
{
}

void EdgeRatioMetric::calculateQuality()
{
    // get all elements of mesh
    const std::vector<MeshLib::Element*>& elements(_mesh.getElements());
    const std::size_t nElements (_mesh.getNumberOfElements());
    for (std::size_t k(0); k < nElements; k++)
    {
        Element const& elem (*elements[k]);
        switch (elem.getGeomType())
        {
        case MeshElemType::LINE:
            _element_quality_metric[k] = 1.0;
            break;
        case MeshElemType::TRIANGLE: {
            _element_quality_metric[k] = checkTriangle(*elem.getNode(0), *elem.getNode(1), *elem.getNode(2));
            break;
        }
        case MeshElemType::QUAD: {
            _element_quality_metric[k] = checkQuad(*elem.getNode(0), *elem.getNode(1), *elem.getNode(2), *elem.getNode(3));
            break;
        }
        case MeshElemType::TETRAHEDRON: {
            _element_quality_metric[k] = checkTetrahedron(*elem.getNode(0), *elem.getNode(1), *elem.getNode(2), *elem.getNode(3));
            break;
        }
        case MeshElemType::PRISM: {
            std::vector<const MathLib::Point3d*> pnts;
            for (std::size_t j(0); j < 6; j++)
                pnts.push_back(elem.getNode(j));
            _element_quality_metric[k] = checkPrism(pnts);
            break;
        }
        case MeshElemType::PYRAMID: {
            std::vector<const MathLib::Point3d*> pnts;
            for (std::size_t j(0); j < 5; j++)
                pnts.push_back(elem.getNode(j));
            _element_quality_metric[k] = checkPyramid(pnts);
            break;
        }
        case MeshElemType::HEXAHEDRON: {
            std::vector<const MathLib::Point3d*> pnts;
            for (std::size_t j(0); j < 8; j++)
                pnts.push_back(elem.getNode(j));
            _element_quality_metric[k] = checkHexahedron(pnts);
            break;
        }
        default:
            ERR ("MeshQualityShortestLongestRatio::check () check for element type %s not implemented.",
                 MeshElemType2String(elem.getGeomType()).c_str());
        }
    }
}

double EdgeRatioMetric::checkTriangle (MathLib::Point3d const& a,
                                       MathLib::Point3d const& b,
                                       MathLib::Point3d const& c) const
{
    double len0(sqrt(MathLib::sqrDist(b, a)));
    double len1(sqrt(MathLib::sqrDist(b, c)));
    double len2(sqrt(MathLib::sqrDist(a, c)));

    if (len0 < len1 && len0 < len2)
    {
        if (len1 < len2)
            return len0 / len2;

        return len0 / len1;
    }

    if (len1 < len2)
    {
        if (len0 < len2)
            return len1 / len2;

        return len1 / len0;
    }

    if (len0 < len1)
        return len2 / len1;

    return len2 / len0;
}

double EdgeRatioMetric::checkQuad (MathLib::Point3d const& a,
                                   MathLib::Point3d const& b,
                                   MathLib::Point3d const& c,
                                   MathLib::Point3d const& d) const
{
    double sqr_lengths[4] = {MathLib::sqrDist (b,a),
                         MathLib::sqrDist (c,b),
                         MathLib::sqrDist (d,c),
                         MathLib::sqrDist (a,d)};

    // sort lengths - since this is a very small array we use bubble sort
    for (std::size_t i(0); i < 4; i++)
        for (std::size_t j(i + 1); j < 4; j++)
            if (sqr_lengths[i] >= sqr_lengths[j])
                std::swap (sqr_lengths[i], sqr_lengths[j]);

    return sqrt(sqr_lengths[0]) / sqrt(sqr_lengths[3]);
}

double EdgeRatioMetric::checkTetrahedron (MathLib::Point3d const& a,
                                          MathLib::Point3d const& b,
                                          MathLib::Point3d const& c,
                                          MathLib::Point3d const& d) const
{
    double sqr_lengths[6] = {MathLib::sqrDist (b,a), MathLib::sqrDist (c,b),
                         MathLib::sqrDist (c,a), MathLib::sqrDist (a,d),
                         MathLib::sqrDist (b,d), MathLib::sqrDist (c,d)};

    // sort lengths - since this is a very small array we use bubble sort
    for (std::size_t i(0); i < 6; i++)
        for (std::size_t j(i + 1); j < 6; j++)
            if (sqr_lengths[i] >= sqr_lengths[j])
                std::swap (sqr_lengths[i], sqr_lengths[j]);

    return sqrt(sqr_lengths[0]) / sqrt(sqr_lengths[5]);
}

double EdgeRatioMetric::checkPrism (std::vector<const MathLib::Point3d*> const& pnts) const
{
    double sqr_lengths[9] = {MathLib::sqrDist (*pnts[0],*pnts[1]),
                         MathLib::sqrDist (*pnts[1],*pnts[2]),
                         MathLib::sqrDist (*pnts[2],*pnts[0]),
                         MathLib::sqrDist (*pnts[3],*pnts[4]),
                         MathLib::sqrDist (*pnts[4],*pnts[5]),
                         MathLib::sqrDist (*pnts[5],*pnts[3]),
                         MathLib::sqrDist (*pnts[0],*pnts[3]),
                         MathLib::sqrDist (*pnts[1],*pnts[4]),
                         MathLib::sqrDist (*pnts[2],*pnts[5])};

    // sort lengths - since this is a very small array we use bubble sort
    for (std::size_t i(0); i < 9; i++)
        for (std::size_t j(i + 1); j < 9; j++)
            if (sqr_lengths[i] >= sqr_lengths[j])
                std::swap (sqr_lengths[i], sqr_lengths[j]);

    return sqrt(sqr_lengths[0]) / sqrt(sqr_lengths[8]);
}

double EdgeRatioMetric::checkPyramid (std::vector<const MathLib::Point3d*> const &pnts) const
{
    double sqr_lengths[8] = {MathLib::sqrDist (*pnts[0],*pnts[1]),
                         MathLib::sqrDist (*pnts[1],*pnts[2]),
                         MathLib::sqrDist (*pnts[2],*pnts[3]),
                         MathLib::sqrDist (*pnts[3],*pnts[0]),
                         MathLib::sqrDist (*pnts[0],*pnts[4]),
                         MathLib::sqrDist (*pnts[1],*pnts[4]),
                         MathLib::sqrDist (*pnts[2],*pnts[4]),
                         MathLib::sqrDist (*pnts[3],*pnts[4])};

    // sort lengths - since this is a very small array we use bubble sort
    for (std::size_t i(0); i < 8; i++)
        for (std::size_t j(i + 1); j < 8; j++)
            if (sqr_lengths[i] >= sqr_lengths[j])
                std::swap (sqr_lengths[i], sqr_lengths[j]);

    return sqrt(sqr_lengths[0]) / sqrt(sqr_lengths[7]);
}

double EdgeRatioMetric::checkHexahedron (std::vector<const MathLib::Point3d*> const & pnts) const
{
    double sqr_lengths[12] = {MathLib::sqrDist (*pnts[0],*pnts[1]),
                          MathLib::sqrDist (*pnts[1],*pnts[2]),
                          MathLib::sqrDist (*pnts[2],*pnts[3]),
                          MathLib::sqrDist (*pnts[3],*pnts[0]),
                          MathLib::sqrDist (*pnts[4],*pnts[5]),
                          MathLib::sqrDist (*pnts[5],*pnts[6]),
                          MathLib::sqrDist (*pnts[6],*pnts[7]),
                          MathLib::sqrDist (*pnts[7],*pnts[4]),
                          MathLib::sqrDist (*pnts[0],*pnts[4]),
                          MathLib::sqrDist (*pnts[1],*pnts[5]),
                          MathLib::sqrDist (*pnts[2],*pnts[6]),
                          MathLib::sqrDist (*pnts[3],*pnts[7])};

    // sort lengths - since this is a very small array we use bubble sort
    for (std::size_t i(0); i < 12; i++)
        for (std::size_t j(i + 1); j < 12; j++)
            if (sqr_lengths[i] >= sqr_lengths[j])
                std::swap (sqr_lengths[i], sqr_lengths[j]);

    return sqrt(sqr_lengths[0]) / sqrt(sqr_lengths[11]);
}
} // end namespace MeshLib
