/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file   RelativePermeabilityCurve.h
 *
 * Created on November 2, 2016, 1:41 PM
 */

#ifndef OGS_RELATIVE_PERMEABILITY_CURVE_H
#define OGS_RELATIVE_PERMEABILITY_CURVE_H

#include <memory>
#include "RelativePermeability.h"

#include "MathLib/MathTools.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"

namespace MaterialLib
{
namespace PorousMedium
{
class RelativePermeabilityCurve final : public RelativePermeability
{
public:
    RelativePermeabilityCurve(
        std::unique_ptr<MathLib::PiecewiseLinearInterpolation>& curve_data)
        : _Sr(curve_data->getSupportMin()),
          _Smax(curve_data->getSupportMax()),
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
        const double S = MathLib::limitValueInInterval(
            saturation, _Sr + _minor_offset, _Smax - _minor_offset);

        return _curve_data->getValue(S);
    }

private:
    const double _Sr;    ///< Residual saturation.
    const double _Smax;  ///< Maximum saturation.

    std::unique_ptr<MathLib::PiecewiseLinearInterpolation> _curve_data;
};
}  // end namespace
}  // end namespace
#endif /* OGS_RELATIVE_PERMEABILITY_CURVE_H */
