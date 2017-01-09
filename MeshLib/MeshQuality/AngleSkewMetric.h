/**
 * \file   AngleSkewMetric.h
 * \author Thomas Fischer
 * \date   2011-03-17
 * \brief  Definition of the AngleSkewMetric class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef ANGLESKEWMETRIC_H_
#define ANGLESKEWMETRIC_H_

#include "ElementQualityMetric.h"

namespace MeshLib
{

/**
 * Calculates the quality of mesh elements based on the EquiAngleSkew measure
 */
class AngleSkewMetric final : public ElementQualityMetric
{
public:
    AngleSkewMetric(Mesh const& mesh);

    void calculateQuality();

private:
    double checkTriangle(Element const& elem) const;
    double checkQuad(Element const& elem) const;
    double checkTetrahedron(Element const& elem) const;
    double checkHexahedron(Element const& elem) const;
    double checkPrism (Element const& elem) const;
    void getMinMaxAngleFromQuad(double const* const n0,
                                double const* const n1, double const* const n2,
                                double const* const n3, double &min_angle,
                                double &max_angle) const;
    void getMinMaxAngleFromTriangle(double const* const n0,
                                    double const* const n1, double const* const n2,
                                    double &min_angle, double &max_angle) const;
};
}

#endif /* ANGLESKEWMETRIC_H_ */
