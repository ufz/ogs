/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on March 4, 2021, 3:05 PM
 */

#pragma once

#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
class Phase;

/**
 * \brief A class for viscosity model that is defined by
 *        The International Association for the Properties of Water and Steam
 *        <a href="http://www.iapws.org/relguide/visc.pdf">IAPWS</a>
 *
 *        With the definition, the viscosity is a function of temperature and
 *        water density
 *
 *  \attention The critical enhancement, \f$\bar{\mu}_2\f$, which is significant
 *             only within the boundaries specified by
 *                 \f[ T (\mbox{in K}) \in (645.91, 650.77) \f]
 *             and
 *                 \f[ \rho (\mbox{in kg m}^{-3}) \in (245.8, 405.3)\f],
 *             is not considered.
 */
class WaterViscosityIAPWS final : public Property
{
public:
    explicit WaterViscosityIAPWS(std::string name) { name_ = std::move(name); }
    void checkScale() const override
    {
        if (!std::holds_alternative<Phase*>(scale_))
        {
            OGS_FATAL(
                "The property 'WaterViscosityIAPWS' is "
                "implemented on the 'Phase' scale only.");
        }
    }

    /// \return The water viscosity.
    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t, double const dt) const override;
    /// \return The derivative of water viscosity with respect to
    /// temperature or phase (water) pressure.
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const variable,
                            ParameterLib::SpatialPosition const& pos,
                            double const t, double const dt) const override;

private:
    static constexpr double ref_T_ = 647.096;  ///< reference temperature in K
    static constexpr double ref_rho_ =
        322.0;  ///< reference density in `kg/m^3`
    static constexpr double ref_mu_ = 1.0e-6;  ///< reference viscosity in Pa.s

    // Coefficients Hi and Hij are given in two static arrays in the cpp file.
};
}  // namespace MaterialPropertyLib
