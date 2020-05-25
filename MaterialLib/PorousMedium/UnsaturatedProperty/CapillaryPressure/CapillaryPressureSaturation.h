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

#include <limits>
#include <string>

namespace MaterialLib
{
namespace PorousMedium
{
/// Base class of capillary pressure models
class CapillaryPressureSaturation
{
public:
    /**
     * @param Sr     Residual saturation, \f$ S_r \f$
     * @param Sg_r     Residual saturation of gas phase, \f$ Sg_r \f$
     * @param Smax   Maximum saturation, \f$ S_{\mbox{max}} \f$
     * @param Pc_max Maximum capillary pressure, \f$ P_c^{\mbox{max}}\f$
     */
    CapillaryPressureSaturation(const double Sr, const double Sg_r,
                                const double Smax, const double Pc_max)
        : saturation_r_(Sr),
          saturation_nonwet_r_(Sg_r),
          saturation_max_(Smax),
          pc_max_(Pc_max)
    {
    }
    virtual ~CapillaryPressureSaturation() = default;

    /// Get model name.
    virtual std::string getName() const = 0;

    /// Get capillary pressure.
    virtual double getCapillaryPressure(const double saturation) const = 0;
    /// Get saturation.
    virtual double getSaturation(const double capillary_pressure) const = 0;

    /// Get the derivative of the capillary pressure with respect to saturation
    virtual double getdPcdS(const double saturation) const = 0;

    /// Get the second derivative of the capillary pressure with respect to
    /// saturation
    virtual double getd2PcdS2(const double saturation) const = 0;

protected:
    const double saturation_r_;         ///< Residual saturation.
    const double saturation_nonwet_r_;  ///< Residual saturation of nonwetting
                                        ///phase (optional).
    const double saturation_max_;       ///< Maximum saturation.
    const double pc_max_;               ///< Maximum capillary pressure

    /** A small number for an offset:
     *  1. to set the bound of S, the saturation, such that
     *     S in  [saturation_r_+minor_offset_, saturation_max_-minor_offset_]
     *  2. to set the bound of Pc, the capillary pressure, such that
     *     Pc in [minor_offset_, pc_max_]
     */
    const double minor_offset_ = std::numeric_limits<double>::epsilon();
};

}  // namespace PorousMedium
}  // namespace MaterialLib
