// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>
#include <string_view>
#include <vector>

namespace MathLib
{
/// Tag selecting the piecewise linear interpolation of a Curve.
struct PiecewiseLinear
{
};

/// A one-dimensional curve, e.g. read from the \c curves section of the
/// project file.
///
/// The interpolation algorithm is hidden behind a pointer to implementation,
/// so users of the class do not depend on it and further interpolation
/// schemes can be added without changing them.
class Curve final
{
public:
    /// Creates a curve interpolating the given supporting points piecewise
    /// linearly.
    /// \param supporting_points  the supporting points \f$x_j\f$, pairwise
    ///     different, not necessarily sorted; duplicates, an empty vector or a
    ///     size mismatch with \c values are a fatal error
    /// \param values  the values \f$y_j\f$ at the supporting points
    Curve(PiecewiseLinear /*tag*/, std::vector<double> supporting_points,
          std::vector<double> values);

    /// As above; the \c curve_name is used in error messages.
    Curve(PiecewiseLinear /*tag*/, std::string_view curve_name,
          std::vector<double> supporting_points, std::vector<double> values);

    Curve(Curve const& other);
    Curve(Curve&& other) noexcept;
    Curve& operator=(Curve const& other);
    Curve& operator=(Curve&& other) noexcept;
    ~Curve();

    /// The value at the given point. Outside of the range of the supporting
    /// points it is constant.
    double getValue(double const x) const;

    /// The derivative at the given point. Outside of the range of the
    /// supporting points it is zero.
    double getDerivative(double const x) const;

    double getSupportMin() const;

    double getSupportMax() const;

private:
    struct Impl;

    std::unique_ptr<Impl> impl_;
};
}  // namespace MathLib
