/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#pragma once

#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
class Medium;
class Phase;
class Component;
/**
 * A simple relative permeability function proposed by
 * Kent S GeneralizedPower \cite udell1985heat.
 *
 *  Definition:
 *          \f[ k_{\text{rel}}
 * =\text{multiplier} \left(1-S^{\text{eff}}\right)^{\text{exponent}}\f] where
 *           \f$S^{\text{eff}}\f$ is the effective saturation of
 * the non-wetting phase.
 *
 * This class handles the non-wetting (gas) phase portion of this relative
 * permeability property for the  relative permeability
 * function.
 *
 * \details This property must be a medium property, it
 * computes the permeability reduction due to saturation as function of
 * capillary pressure.
 */
class RelPermGeneralizedPowerNonwettingPhase final : public Property
{
private:
    const double residual_liquid_saturation_;
    const double residual_gas_saturation_;
    const double min_relative_permeability_;
    const double a_;
    const double lambda_;

public:
    RelPermGeneralizedPowerNonwettingPhase(
        std::string name,
        const double residual_liquid_saturation,
        const double residual_gas_saturation,
        const double min_relative_permeability,
        const double a,
        const double lambda);

    void checkScale() const override
    {
        if (!std::holds_alternative<Medium*>(scale_))
        {
            OGS_FATAL(
                "The property "
                "'RelativePermeabilityGeneralizedPowerNonwettingPhase' is "
                "implemented on the 'media' scale only.");
        }
    }

    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t, double const dt) const override;
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const variable,
                            ParameterLib::SpatialPosition const& pos,
                            double const t, double const dt) const override;
};

}  // namespace MaterialPropertyLib
