/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   CapillaryPressureSaturationCurve.h
 *
 * Created on November 3, 2016, 9:50 AM
 */

#pragma once

#include <memory>

#include "CapillaryPressureSaturation.h"

#include "MathLib/Curve/PiecewiseLinearMonotonicCurve.h"
#include "MathLib/MathTools.h"

namespace MaterialLib
{
namespace PorousMedium
{
class CapillaryPressureSaturationCurve final
    : public CapillaryPressureSaturation
{
public:
    CapillaryPressureSaturationCurve(
        std::unique_ptr<MathLib::PiecewiseLinearMonotonicCurve>&& curve_data)
        : CapillaryPressureSaturation(
              curve_data->getSupportMin(),
              1.0 - curve_data->getSupportMax(),
              curve_data->getSupportMax(),
              curve_data->getValue(curve_data->getSupportMin())),
          _curve_data(std::move(curve_data))
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
        const double S = MathLib::limitValueInInterval(
            saturation, _saturation_r + _minor_offset,
            _saturation_max - _minor_offset);

        return _curve_data->getValue(S);
    }

    /// Get saturation.
    double getSaturation(const double capillary_pressure) const override
    {
        const double pc = MathLib::limitValueInInterval(capillary_pressure,
                                                        _minor_offset, _pc_max);
        return _curve_data->getInverseVariable(pc);
    }

    /// Get the derivative of the capillary pressure with respect to saturation
    double getdPcdS(const double saturation) const override
    {
        const double S = MathLib::limitValueInInterval(
            saturation, _saturation_r + _minor_offset,
            _saturation_max - _minor_offset);

        return _curve_data->getDerivative(S);
    }

    /// Get the second derivative of the capillary pressure with respect to
    /// saturation
    /// In the case of piecewise linear curves, it is always zero.
    double getd2PcdS2(const double /*saturation*/) const override { return 0; }

private:
    std::unique_ptr<MathLib::PiecewiseLinearMonotonicCurve> _curve_data;
};

}  // end namespace
}  // end namespace
