/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file:   RelativePermeability.h
 *
 */

#ifndef OGS_RELATIVE_PERMEABILITY_H
#define OGS_RELATIVE_PERMEABILITY_H

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
    /** A small number for an offset:
     *  1. to set the bound of S, the saturation, such that
     *     S in  [_Sr+_minor_offset, _Smax-_minor_offset]
     *  2. to set the bound of Pc, the capillary pressure, such that
     *     Pc in [_minor_offset, _Pc_max]
     */
    const double _minor_offset = std::numeric_limits<double>::epsilon();
};

}  // end namespace
}  // end namespace

#endif /* OGS_RELATIVE_PERMEABILITY_H */
