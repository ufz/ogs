/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file:   RelativePermeability.h
 *
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
        : _saturation_r(Sr), _saturation_max(Smax)
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
    /// that S in  [Sr+_minor_offset, Smax-_minor_offset].
    const double _minor_offset = std::numeric_limits<double>::epsilon();

    const double _saturation_r;    ///< Residual saturation.
    const double _saturation_max;  ///< Maximum saturation.
};

}  // end namespace
}  // end namespace
