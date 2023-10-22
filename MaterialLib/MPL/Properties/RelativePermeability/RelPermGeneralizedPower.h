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
 * A generalized power law for the relative permeability
 *
 *  Definition:
 *          \f[ k_{\text{rel}}
 * = \text{multiplier} {S_{\text{eff}}}^{\text{exponent}}\f] where
 *          \f$S_{\text{eff}}\f$ is the effective saturation of
 * liquid phase.
 *
 * This class handles the wetting (liquid) phase portion of this relative
 * permeability property.
 *
 * \details This property must be a medium property, it
 * computes the permeability reduction due to saturation as function of
 * capillary pressure.
 */
class RelPermGeneralizedPower final : public Property
{
private:
    const double residual_liquid_saturation_;
    const double residual_gas_saturation_;
    const double min_relative_permeability_;
    const double a_;
    const double lambda_;

public:
    RelPermGeneralizedPower(std::string name,
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
                "The property 'RelativePermeabilityGeneralizedPower' is "
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
