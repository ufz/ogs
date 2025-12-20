// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
