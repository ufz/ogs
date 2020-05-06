/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#pragma once

#include "BaseLib/ConfigTree.h"
#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
class Medium;
class Phase;
class Component;

/// Van Genuchten relative permeability function.
///
/// This property must be a medium property, it computes the permeability
/// reduction due to saturation as function of capillary pressure.
class RelPermVanGenuchten final : public Property
{
private:
    double const S_L_res_;
    double const S_L_max_;
    double const k_rel_min_;
    double const m_;
    Medium* medium_ = nullptr;

public:
    RelPermVanGenuchten(double const residual_liquid_saturation,
                        double const residual_gas_saturation,
                        double const min_relative_permeability_liquid,
                        double const exponent);
    /// This method assigns a pointer to the material object that is the owner
    /// of this property
    void setScale(
        std::variant<Medium*, Phase*, Component*> scale_pointer) override
    {
        if (std::holds_alternative<Medium*>(scale_pointer))
        {
            medium_ = std::get<Medium*>(scale_pointer);
        }
        else
        {
            OGS_FATAL(
                "The property 'RelativePermeabilityVanGenuchten' is "
                "implemented on the 'media' scale only.");
        }
    }

    /// Those methods override the base class implementations and
    /// actually compute and set the property values_ and dValues_.
    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t, double const dt) const override;
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const variable,
                            ParameterLib::SpatialPosition const& pos,
                            double const t, double const dt) const override;
};
}  // namespace MaterialPropertyLib
