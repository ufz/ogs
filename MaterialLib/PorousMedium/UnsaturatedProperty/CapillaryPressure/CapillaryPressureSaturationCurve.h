/**
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file
 *
 * Created on November 3, 2016, 9:50 AM
 */

#pragma once

#include <algorithm>
#include <memory>

#include "CapillaryPressureSaturation.h"
#include "MathLib/Curve/PiecewiseLinearMonotonicCurve.h"

namespace MaterialLib
{
namespace PorousMedium
{
class CapillaryPressureSaturationCurve final
    : public CapillaryPressureSaturation
{
public:
    explicit CapillaryPressureSaturationCurve(
        std::unique_ptr<MathLib::PiecewiseLinearMonotonicCurve>&& curve_data)
        : CapillaryPressureSaturation(
              curve_data->getSupportMin(),
              1.0 - curve_data->getSupportMax(),
              curve_data->getSupportMax(),
              curve_data->getValue(curve_data->getSupportMin())),
          curve_data_(std::move(curve_data))
    {
    }

    /// Get model name.
    std::string getName() const override
    {
        return "Capillary pressure saturation curve.";
    }

    /// Get capillary pressure.
    double getCapillaryPressure(const double saturation) const override
    {
        const double S = std::clamp(saturation, saturation_r_ + minor_offset_,
                                    saturation_max_ - minor_offset_);

        return curve_data_->getValue(S);
    }

    /// Get saturation.
    double getSaturation(const double capillary_pressure) const override
    {
        const double pc =
            std::clamp(capillary_pressure, minor_offset_, pc_max_);
        return curve_data_->getInverseVariable(pc);
    }

    /// Get the derivative of the capillary pressure with respect to saturation
    double getdPcdS(const double saturation) const override
    {
        const double S = std::clamp(saturation, saturation_r_ + minor_offset_,
                                    saturation_max_ - minor_offset_);

        return curve_data_->getDerivative(S);
    }

    /// Get the second derivative of the capillary pressure with respect to
    /// saturation
    /// In the case of piecewise linear curves, it is always zero.
    double getd2PcdS2(const double /*saturation*/) const override { return 0; }

private:
    std::unique_ptr<MathLib::PiecewiseLinearMonotonicCurve> curve_data_;
};

}  // namespace PorousMedium
}  // namespace MaterialLib
