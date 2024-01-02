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
 * \brief A class for thermal conductivity model that is defined by
 *        The International Association for the Properties of Water and Steam
 *        <a href="http://www.iapws.org/relguide/ThCond.pdf">IAPWS</a>
 *        (File accessed at 13.01.2023) - (Daucik and Dooley, 2011)
 *
 *        With the definition, the thermal conductivity is a function of
 * temperature and water density
 *
 *  \attention The critical enhancement, \f$\bar{\lambda}_2\f$, is not
 * considered. For information on region of significance and the significance,
 *              please see Figure 2 from the document linked in the upper
 * paragraph.
 */
class WaterThermalConductivityIAPWS final : public Property
{
public:
    explicit WaterThermalConductivityIAPWS(std::string name)
    {
        name_ = std::move(name);
    }
    void checkScale() const override
    {
        if (!std::holds_alternative<Phase*>(scale_))
        {
            OGS_FATAL(
                "The property 'WaterThermalConductivityIAPWS' is "
                "implemented on the 'Phase' scale only.");
        }
    }

    /// \return The water thermal conductivity.
    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t, double const dt) const override;
    /// \return The derivative of water thermal conductivity with respect to
    /// temperature or phase (water) pressure.
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const variable,
                            ParameterLib::SpatialPosition const& pos,
                            double const t, double const dt) const override;

private:
    static constexpr double ref_T_ = 647.096;  ///< reference temperature in K
    static constexpr double ref_rho_ =
        322.0;  ///< reference density in `kg/m^3`
    static constexpr double ref_lambda_ =
        1.0e-3;  ///< reference thermal conductivity in `W.K^-1.m^-1`
};
}  // namespace MaterialPropertyLib
