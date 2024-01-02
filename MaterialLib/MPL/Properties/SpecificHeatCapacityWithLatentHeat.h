/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on May 20, 2022
 */
#pragma once

#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
/**
 * Model for the apparent specific heat capacity including latent heat
 *\f$\ell\f$
 *
 * \details This model is for media with a phase change. This property must be
 * a medium property, it computes the apparent specific heat capacity based on
 * a phase transition spread over a temperature interval:
 *
 * \f[
 *      C_{\mathrm{app}} = \left[C - \varrho_\mathrm{fR} \ell \frac{\partial
 *                         \phi_\mathrm{f}}{\partial T} \right]
 * \quad\rightarrow\quad c_{\mathrm{app}} = \frac{C_{\mathrm{app}}}{\varrho}
 * \f]
 *
 * with \f$\ell\f$ as specific enthalpy of melting, \f$C\f$ as effective
 * volumetric heat capacity and \f$\varrho\f$ as the effective density of the
 * mixture, \f$\varrho_\mathrm{fR}\f$ as the real density of the frozen phase
 * as well as \f$\phi_\mathrm{f}\f$ as the temperature-dependent frozen
 * liquid volume fraction.
 **/
class SpecificHeatCapacityWithLatentHeat final : public Property
{
    struct PhaseProperties
    {
        Property const* liquid = nullptr;
        Property const* frozen = nullptr;
        Property const* porous = nullptr;
    };

public:
    SpecificHeatCapacityWithLatentHeat(std::string name, double const l);

    void checkScale() const override;

    // initialize container with pointers to the phases properties
    void setProperties(
        std::vector<std::unique_ptr<Phase>> const& phases) override;

    double effectiveVolumetricHeatCapacity(
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
    double const l_;  //< specific latent heat of melting

    /// Pointers to the properties in each phase.
    PhaseProperties densities_;
    /// Pointers to the properties in each phase.
    PhaseProperties spec_heat_capacities_;
};
}  // namespace MaterialPropertyLib
