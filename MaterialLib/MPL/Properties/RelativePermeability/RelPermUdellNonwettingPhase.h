/**
 * \file
 * \author Norbert Grunwald
 * \date   01.12.2020
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
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
 * Kent S Udell \cite udell1985heat.
 *
 *  Definition:
 *          \f[ k_{\text{rel}}^{\alpha}
 * =\left(S^{\text{eff}}_{\alpha}\right)^{3}\f] where
 *          - \f$k_{\text{rel}}^{\alpha}\f$ is relative permeability of phase
 * \f$\alpha\f$
 *          - \f$S^{\text{eff}}_{\alpha}\f$ is the effective saturation of
 * phase \f$\alpha\f$
 *
 * This class handles the non-wetting (gas) phase portion of this relative
 * permeability property, i,e with \f$\alpha=g\f$ for the  relative permeability
 * function.
 *
 * \details This property must be a medium property, it
 * computes the permeability reduction due to saturation as function of
 * capillary pressure.
 */
class RelPermUdellNonwettingPhase final : public Property
{
private:
    const double residual_liquid_saturation_;
    const double residual_gas_saturation_;
    const double min_relative_permeability_;

public:
    RelPermUdellNonwettingPhase(std::string name,
                                const double residual_liquid_saturation,
                                const double residual_gas_saturation,
                                const double min_relative_permeability);

    void checkScale() const override
    {
        if (!std::holds_alternative<Medium*>(scale_))
        {
            OGS_FATAL(
                "The property 'RelativePermeabilityUdellNonwettingPhase' is "
                "implemented on the 'media' scale only.");
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
