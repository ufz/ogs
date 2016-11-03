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
    /// \param saturation Non-wetting phase saturation
    virtual double getValue(const double saturation) const = 0;
};

}  // end namespace
}  // end namespace

#endif /* OGS_RELATIVE_PERMEABILITY_H */
