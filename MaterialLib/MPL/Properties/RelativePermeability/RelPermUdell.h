/**
 * \file
 * \author Norbert Grunwald
 * \date   01.12.2020
 * \brief
 *
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
/**
 * \class RelPermUdell
 * \brief Relative permeability used in:
 *        Kent S Udell. "Heat transfer in porous media considering phase change
 *        and capillarity-the heat pipe effect".In:International Journal of Heat
 *        and Mass Transfer 28.2 (1985), pp. 485-495
 * \details This property must be a medium
 * property, it computes the permeability reduction due to saturation as
 * function of capillary pressure.
 */
class RelPermUdell final : public Property
{
private:
    const double residual_liquid_saturation_;
    const double residual_gas_saturation_;
    const double min_relative_permeability_liquid_;
    const double min_relative_permeability_gas_;

public:
    RelPermUdell(std::string name, const double /*residual_liquid_saturation*/,
                 const double /*residual_gas_saturation*/,
                 const double /*min_relative_permeability_liquid_*/,
                 const double /*min_relative_permeability_gas_*/);

    void checkScale() const override
    {
        if (!std::holds_alternative<Medium*>(scale_))
        {
            OGS_FATAL(
                "The property 'RelPermUdell' is implemented on the "
                "'media' scale only.");
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
