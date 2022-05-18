/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */
#pragma once

#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
/**
 * Model for the apparent specific heat capacity including latent heat $\ell$
 *
 * \details This model is for media with a phase change. This property must be
 * a medium property, it computes the apparent specific heat capacity based on
 * a phase transition spread over a temperature interval.
 *
 * \f$
 *      C_{\mathrm{app}} = \left[C - \varrho \ell \frac{\partial
 *                         \phi_\text{f}}{\partial T} \right]
 * \f$
 *
 * with $\ell$ as specific enthalpy of melting, $C$ as effective volumetric heat
 * capacity and $\varrho$ as the effective density of the mixture, as well as
 * $\phi_\text{f}$ as the temperature-dependant frozen liquid volume fraction.
 **/
class SpecificHeatCapacityWithLatentHeat final : public Property
{
    struct PhaseProperties
    {
        const Property* liquid;
        const Property* frozen;
        const Property* porous;
    };

public:
    SpecificHeatCapacityWithLatentHeat(std::string name, double const L);

    void checkScale() const override;

    // initialize container with pointers to the phases properties
    void setProperties(
        std::vector<std::unique_ptr<Phase>> const& phases) override;

    double mixtureVolumetricHeatCapacity(
        VariableArray const& variable_array,
        ParameterLib::SpatialPosition const& pos,
        double const t,
        double const dt) const;

    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t,
                           double const dt) const override;

    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const variable,
                            ParameterLib::SpatialPosition const& pos,
                            double const t,
                            double const dt) const override;

private:
    double const L_;  //< Latent heat of melting (enthalpy of the first order
                      // phase change).

    /// Pointers to the properties in each phase.
    PhaseProperties densities_;
    /// Pointers to the properties in each phase.
    PhaseProperties spec_heat_capacities_;
};
}  // namespace MaterialPropertyLib
