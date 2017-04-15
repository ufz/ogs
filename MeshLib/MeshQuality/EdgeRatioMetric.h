/**
 * \file   EdgeRatioMetric.h
 * \author Thomas Fischer
 * \date   2011-03-03
 * \brief  Definition of the AreaMetric class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ElementQualityMetric.h"
#include "MathLib/Point3d.h"

namespace MeshLib
{

/**
 * Calculates the quality of mesh elements based on the ratio between shortest and longest edge of an element
 */
class EdgeRatioMetric : public ElementQualityMetric
{
public:
    EdgeRatioMetric(Mesh const& mesh);
    virtual ~EdgeRatioMetric() = default;

    virtual void calculateQuality ();

private:
    double checkTriangle (MathLib::Point3d const& a,
                          MathLib::Point3d const& b,
                          MathLib::Point3d const& c) const;
    double checkQuad (MathLib::Point3d const& a,
                      MathLib::Point3d const& b,
                      MathLib::Point3d const& c,
                      MathLib::Point3d const& d) const;
    double checkTetrahedron (MathLib::Point3d const& a,
                             MathLib::Point3d const& b,
                             MathLib::Point3d const& c,
                             MathLib::Point3d const& d) const;
    double checkPrism (std::vector<const MathLib::Point3d*> const& pnts) const;
    double checkPyramid (std::vector<const MathLib::Point3d*> const& pnts) const;
    double checkHexahedron (std::vector<const MathLib::Point3d*> const& pnts) const;
};
}
