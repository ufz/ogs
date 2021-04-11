/**
 * \file
 * \author Thomas Fischer
 * \date   2010-01-12
 * \brief  Definition of the Point class.
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

// GeoLib
#include "GeoObject.h"

// MathLib
#include "MathLib/Point3dWithID.h"

namespace GeoLib
{
/**
 * \ingroup GeoLib
 */

// forward declaration
class PointVec;

class Point : public MathLib::Point3dWithID, public GeoLib::GeoObject
{
public:
    Point() = default;

    Point(double x1, double x2, double x3,
          std::size_t id = std::numeric_limits<std::size_t>::max())
        : MathLib::Point3dWithID(std::array<double, 3>({{x1, x2, x3}}), id)
    {
    }

    Point(MathLib::Point3d const& x, std::size_t id)
        : MathLib::Point3dWithID(x, id)
    {
    }

    explicit Point(std::array<double, 3> const& x,
                   std::size_t id = std::numeric_limits<std::size_t>::max())
        : MathLib::Point3dWithID(x, id)
    {
    }

    /// return a geometry type
    GEOTYPE getGeoType() const override { return GEOTYPE::POINT; }

protected:
    friend PointVec;
};

}  // namespace GeoLib
