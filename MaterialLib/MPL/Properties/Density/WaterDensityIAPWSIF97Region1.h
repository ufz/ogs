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

#include "MaterialLib/Fluid/GibbsFreeEnergy/DimensionLessGibbsFreeEnergyRegion1.h"
#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
class Phase;

/** Water density model
 *  base on the IAPWS Industrial Formulation 1997
 *  <a href="http://www.iapws.org/relguide/IF97-Rev.pdf">IF97-Rev</a>
 */
class WaterDensityIAPWSIF97Region1 final : public Property
{
public:
    explicit WaterDensityIAPWSIF97Region1(std::string name)
        : gibbs_free_energy_(
              MaterialLib::Fluid::DimensionLessGibbsFreeEnergyRegion1())
    {
        name_ = std::move(name);
    }
    void checkScale() const override
    {
        if (!std::holds_alternative<Phase*>(scale_))
        {
            OGS_FATAL(
                "The property 'WaterDensityIAPWSIF97Region1' is "
                "implemented on the 'Phase' scale only.");
        }
    }

    /// \return The water density.
    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t, double const dt) const override;
    /// \return The derivative of  water density with respect to
    /// temperature or phase (water) pressure.
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const variable,
                            ParameterLib::SpatialPosition const& pos,
                            double const t, double const dt) const override;

private:
    const MaterialLib::Fluid::DimensionLessGibbsFreeEnergyRegion1
        gibbs_free_energy_;

    static constexpr double ref_T_ = 1386;     ///< reference temperature in K.
    static constexpr double ref_p_ = 1.653e7;  ///< reference pressure in Pa.
};
}  // namespace MaterialPropertyLib
