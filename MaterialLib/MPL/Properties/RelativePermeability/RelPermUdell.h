/**
 * \file
 * \author Norbert Grunwald
 * \date   01.12.2020
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
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
 * A simple relative permeability function as used in
 *
 *        Kent S Udell. "Heat transfer in porous media considering phase change
 *        and capillarity-the heat pipe effect".In:International Journal of Heat
 *        and Mass Transfer 28.2 (1985), pp. 485-495
 *         Definition:
 *          \f$k_{\textrm{rel}}^{\alpha} =
 * \left(s^{\textrm{eff}}_{\alpha}\right)^{3}\f$ where
 *          - \f$k_{\textrm{rel}}^{\alpha}\f$ is relative permeability of phase
 * \f$\alpha\f$
 *          - \f$s^{\textrm{eff}}_{\alpha}\f$ is the effective saturation of
 * phase \f$\alpha\f$ \details This property must be a medium property, it
 * computes the permeability reduction due to saturation as function of
 * capillary pressure.
 */
class RelPermUdell final : public Property
{
private:
    const double residual_liquid_saturation_;
    const double residual_gas_saturation_;
    const double min_relative_permeability_liquid_;
    const double min_relative_permeability_gas_;

public:
    RelPermUdell(std::string name, const double residual_liquid_saturation,
                 const double residual_gas_saturation,
                 const double min_relative_permeability_liquid,
                 const double min_relative_permeability_gas);

    void checkScale() const override
    {
        if (!std::holds_alternative<Medium*>(scale_))
        {
            OGS_FATAL(
                "The property 'RelPermUdell' is implemented on the "
                "'media' scale only.");
        }
    }

    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t, double const dt) const override;
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const primary_variable,
                            ParameterLib::SpatialPosition const& pos,
                            double const t, double const dt) const override;
};

}  // namespace MaterialPropertyLib
