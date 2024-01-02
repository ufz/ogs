/**
 * \file
 * \author Norbert Grunwald
 * \date   27.06.2018
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#pragma once

#include <string>

#include "BaseLib/ConfigTree-fwd.h"
#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
class Medium;
class Phase;
class Component;
/**
 * \class RelPermBrooksCoreyNonwettingPhase
 * \brief Relative permeability function of the non-wetting phase proposed by
 * Brooks&Corey.
 *
 * \details This property must be a medium property, it
 * computes the permeability reduction due to saturation as function
 * of capillary pressure.
 */
class RelPermBrooksCoreyNonwettingPhase final : public Property
{
private:
    const double residual_liquid_saturation_;
    const double residual_gas_saturation_;
    const double min_relative_permeability_;
    const double exponent_;

public:
    RelPermBrooksCoreyNonwettingPhase(std::string name,
                                      const double residual_liquid_saturation,
                                      const double residual_gas_saturation,
                                      const double min_relative_permeability,
                                      const double exponent);

    void checkScale() const override
    {
        if (!std::holds_alternative<Medium*>(scale_))
        {
            OGS_FATAL(
                "The property 'RelPermBrooksCoreyNonwettingPhase' is "
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
