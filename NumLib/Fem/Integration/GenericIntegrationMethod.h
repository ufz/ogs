/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <vector>

#include "MathLib/WeightedPoint.h"

namespace NumLib
{
//! Provides data for arbitrary numerical integration methods of any integration
//! order.
//!
//! This class basically holds a collection of integration points and
//! integration weights.
class GenericIntegrationMethod final
{
public:
    GenericIntegrationMethod(unsigned const order,
                             std::vector<MathLib::WeightedPoint>&& points)
        : order_{order}, points_{std::move(points)}
    {
    }

    unsigned getIntegrationOrder() const { return order_; }

    unsigned getNumberOfPoints() const { return points_.size(); }

    MathLib::WeightedPoint const& getWeightedPoint(unsigned const igp) const
    {
        return points_[igp];
    }

    // forbid accidental copies
    GenericIntegrationMethod(GenericIntegrationMethod const&) = delete;
    GenericIntegrationMethod& operator=(GenericIntegrationMethod const&) =
        delete;
    GenericIntegrationMethod(GenericIntegrationMethod&&) = default;
    GenericIntegrationMethod& operator=(GenericIntegrationMethod&&) = default;

private:
    unsigned order_;

    std::vector<MathLib::WeightedPoint> points_;
};
}  // namespace NumLib
