// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "Curve.h"

#include <cassert>
#include <utility>
#include <variant>

#include "BaseLib/Error.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"

namespace MathLib
{
struct Curve::Impl
{
    std::variant<PiecewiseLinearInterpolation> interpolation;
};

Curve::Curve(PiecewiseLinear /*tag*/, std::vector<double> supporting_points,
             std::vector<double> values)
{
    if (supporting_points.empty() || values.empty())
    {
        OGS_FATAL(
            "Creating a curve from an empty list of supporting points or "
            "values.");
    }
    if (supporting_points.size() != values.size())
    {
        OGS_FATAL(
            "Creating a curve from {:d} supporting points and {:d} values; "
            "the two sizes must match.",
            supporting_points.size(), values.size());
    }

    impl_ = std::make_unique<Impl>(PiecewiseLinearInterpolation{
        std::move(supporting_points), std::move(values)});
}

Curve::Curve(Curve const& other)
{
    // A copy is made once per curve, so the check is cheaper than the copy.
    if (!other.impl_)
    {
        OGS_FATAL("Copying a moved-from Curve.");
    }
    impl_ = std::make_unique<Impl>(*other.impl_);
}

Curve::Curve(Curve&& other) noexcept = default;

Curve& Curve::operator=(Curve const& other)
{
    // The new interpolant is built before the old one is released, so
    // assigning a curve to itself is safe.
    if (!other.impl_)
    {
        OGS_FATAL("Assigning from a moved-from Curve.");
    }
    impl_ = std::make_unique<Impl>(*other.impl_);
    return *this;
}

Curve& Curve::operator=(Curve&& other) noexcept = default;

Curve::~Curve() = default;

// The getters may be called very frequently, so they keep the debug-only
// check instead of the formatting and throwing of OGS_FATAL.
double Curve::getValue(double const x) const
{
    assert(impl_);
    return std::visit([x](auto const& i) { return i.getValue(x); },
                      impl_->interpolation);
}

double Curve::getDerivative(double const x) const
{
    assert(impl_);
    return std::visit([x](auto const& i) { return i.getDerivative(x); },
                      impl_->interpolation);
}

double Curve::getSupportMin() const
{
    assert(impl_);
    return std::visit([](auto const& i) { return i.getSupportMin(); },
                      impl_->interpolation);
}

double Curve::getSupportMax() const
{
    assert(impl_);
    return std::visit([](auto const& i) { return i.getSupportMax(); },
                      impl_->interpolation);
}
}  // namespace MathLib
