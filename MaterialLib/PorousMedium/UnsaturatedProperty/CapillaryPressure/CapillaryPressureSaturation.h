/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file:   CapillaryPressureSaturation.h
 *
 */

#ifndef OGS_CAPILLARY_PRESSURE_SATURATION_H
#define OGS_CAPILLARY_PRESSURE_SATURATION_H

#include <string>

namespace MaterialLib
{
namespace PorousMedium
{
/// Base class of capillary pressure models
class CapillaryPressureSaturation
{
public:
    virtual ~CapillaryPressureSaturation() = default;

    /// Get model name.
    virtual std::string getName() const = 0;

    /// Get capillary pressure.
    virtual double getCapillaryPressure(const double saturation) const = 0;

    /// Get capillary pressure.
    virtual double getSturation(const double capillary_ressure) const = 0;

    /// Get the derivative of the capillary pressure with respect to saturation
    virtual double getdPcdS(const double saturation) const = 0;
};

}  // end namespace
}  // end namespace

#endif /* OGS_CAPILLARY_PRESSURE_SATURATION_H */
