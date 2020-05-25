/**
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file
 */

#pragma once

#include <string>
#include <limits>

namespace MaterialLib
{
namespace PorousMedium
{
/// Base class of relative permeability models
class RelativePermeability
{
public:
     /// @param Sr       Residual saturation.
     /// @param Smax     Maximum saturation.
    RelativePermeability(const double Sr, const double Smax)
        : saturation_r_(Sr), saturation_max_(Smax)
    {
    }

    virtual ~RelativePermeability() = default;

    /// Get model name.
    virtual std::string getName() const = 0;

    /// Get relative permeability value.
    /// \param saturation Wetting phase saturation
    virtual double getValue(const double saturation) const = 0;

    /// Get the derivative of relative permeability with respect to saturation.
    /// \param saturation Wetting phase saturation
    virtual double getdValue(const double saturation) const = 0;

protected:
    /// A small number for an offset to set the bound of S, the saturation, such
    /// that S in  [Sr+minor_offset_, Smax-minor_offset_].
    const double minor_offset_ = std::numeric_limits<double>::epsilon();

    const double saturation_r_;    ///< Residual saturation.
    const double saturation_max_;  ///< Maximum saturation.
};

}  // namespace PorousMedium
}  // namespace MaterialLib
