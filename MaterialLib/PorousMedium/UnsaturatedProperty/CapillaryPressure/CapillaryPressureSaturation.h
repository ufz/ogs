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
    /**
     * @param Sr     Residual saturation, \f$ S_r \f$
     * @param Smax   Maximum saturation, \f$ S_{\mbox{max}} \f$
     * @param Pc_max Maximum capillary pressure, \f$ P_c^{\mbox{max}}\f$
     */
    CapillaryPressureSaturation(const double Sr, const double Smax,
                                const double Pc_max)
        : _saturation_r(Sr), _saturation_max(Smax), _pc_max(Pc_max)
    {
    }
    virtual ~CapillaryPressureSaturation() = default;

    /// Get model name.
    virtual std::string getName() const = 0;

    /// Get capillary pressure.
    virtual double getCapillaryPressure(const double saturation) const = 0;

    /// Get saturation.
    virtual double getSaturation(const double capillary_ressure) const = 0;

    /// Get the derivative of the capillary pressure with respect to saturation
    virtual double getdPcdS(const double saturation) const = 0;

protected:
    const double _saturation_r;    ///< Residual saturation.
    const double _saturation_max;  ///< Maximum saturation.
    const double _pc_max;          ///< Maximum capillaray pressure

    /** A small number for an offset:
     *  1. to set the bound of S, the saturation, such that
     *     S in  [_saturation_r+_minor_offset, _saturation_max-_minor_offset]
     *  2. to set the bound of Pc, the capillary pressure, such that
     *     Pc in [_minor_offset, _pc_max]
     */
    const double _minor_offset = std::numeric_limits<double>::epsilon();
};

}  // end namespace
}  // end namespace

#endif /* OGS_CAPILLARY_PRESSURE_SATURATION_H */
