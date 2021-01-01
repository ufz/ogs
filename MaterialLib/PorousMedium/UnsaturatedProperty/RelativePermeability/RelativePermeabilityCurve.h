/**
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file
 *
 * Created on November 2, 2016, 1:41 PM
 */

#pragma once

#include <algorithm>
#include <memory>

#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "RelativePermeability.h"

namespace MaterialLib
{
namespace PorousMedium
{
class RelativePermeabilityCurve final : public RelativePermeability
{
public:
    explicit RelativePermeabilityCurve(
        std::unique_ptr<MathLib::PiecewiseLinearInterpolation>&& curve_data)
        : RelativePermeability(curve_data->getSupportMin(),
                               curve_data->getSupportMax()),
          _curve_data(std::move(curve_data))
    {
    }

    /// Get model name.
    std::string getName() const override
    {
        return "Relative permeability curve.";
    }

    /// Get relative permeability value.
    double getValue(const double saturation) const override
    {
        const double S = std::clamp(saturation, _saturation_r + _minor_offset,
                                    _saturation_max - _minor_offset);

        return _curve_data->getValue(S);
    }

    /// Get the derivative of relative permeability with respect to saturation.
    /// \param saturation Wetting phase saturation
    double getdValue(const double saturation) const override
    {
        const double S = std::clamp(saturation, _saturation_r + _minor_offset,
                                    _saturation_max - _minor_offset);

        return _curve_data->getDerivative(S);
    }

private:
    std::unique_ptr<MathLib::PiecewiseLinearInterpolation> _curve_data;
};
}  // namespace PorousMedium
}  // namespace MaterialLib
