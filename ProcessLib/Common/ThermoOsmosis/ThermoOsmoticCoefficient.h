// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "BaseLib/Error.h"
#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Utils/FormEigenTensor.h"
#include "ParameterLib/SpatialPosition.h"

namespace ProcessLib
{
/// Returns the thermo-osmotic coefficient tensor \f$k_T\f$
/// (\f$[k_T] = m^2/(K \cdot s)\f$) of the given medium.
///
/// The coefficient can be parametrised either directly via the
/// \c thermal_osmosis_coefficient medium property
/// (\f$[k_T] = m^2/(K \cdot s)\f$), or indirectly via the
/// \c thermal_osmosis_permeability medium property, in which case
/// \f$k_T = \epsilon_T k / \mu\f$ with the scalar thermo-osmotic permeability
/// \f$\epsilon_T\f$ (\f$[\epsilon_T] = Pa/K\f$), which converts a temperature
/// gradient into the equivalent pore pressure gradient driving the
/// thermo-osmotic flux, the intrinsic permeability
/// \f$k\f$ (\f$[k] = m^2\f$), and the liquid's dynamic viscosity
/// \f$\mu\f$ (\f$[\mu] = Pa \cdot s\f$).
/// If neither is defined, the zero tensor is returned (no thermo-osmosis).
/// That the two properties are not defined at the same time, and that neither
/// is left on the solid phase, is checked once at process creation by
/// checkThermoOsmosisProperties() and not re-checked here.
///
/// \f$k\f$ and \f$\mu\f$ are passed in rather than read from the medium: the
/// callers have already evaluated them for the Darcy term, and the caller's
/// values are the ones evaluated with a fully populated variable array (in
/// ThermoRichardsMechanics the array reaching this function holds no primary
/// variables at all, so a temperature-dependent \f$\mu\f$ read here would
/// evaluate to NaN). The same limitation applies to the two thermo-osmosis
/// properties themselves, which are read with that array: in
/// ThermoRichardsMechanics they must not depend on primary variables.
///
/// \note \f$\epsilon_T\f$ is a scalar, so \f$k_T\f$ inherits the anisotropy of
/// \f$k\f$ alone and the composition order does not arise. A tensor-valued
/// \c thermal_osmosis_permeability is rejected by the property's value access.
///
/// @tparam GlobalDim spatial dimension of the returned tensor
/// @param medium the medium the thermo-osmosis properties are read from
/// @param variable_array primary variables at the evaluation point
/// @param pos spatial position of the evaluation point
/// @param t current time
/// @param dt current time increment
/// @param intrinsic_permeability the intrinsic permeability tensor \f$k\f$
/// (\f$[k] = m^2\f$) at the evaluation point
/// @param liquid_dynamic_viscosity the liquid's dynamic viscosity \f$\mu\f$
/// (\f$[\mu] = Pa \cdot s\f$) at the evaluation point
template <int GlobalDim>
Eigen::Matrix<double, GlobalDim, GlobalDim> getThermoOsmoticCoefficient(
    MaterialPropertyLib::Medium const& medium,
    MaterialPropertyLib::VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& pos,
    double const t,
    double const dt,
    Eigen::Matrix<double, GlobalDim, GlobalDim> const& intrinsic_permeability,
    double const liquid_dynamic_viscosity)
{
    if (medium.hasProperty(
            MaterialPropertyLib::PropertyType::thermal_osmosis_permeability))
    {
        if (liquid_dynamic_viscosity <= 0.)
        {
            OGS_FATAL(
                "Liquid dynamic viscosity must be > 0 when using "
                "thermal_osmosis_permeability, but is {:g} in {:s}.",
                liquid_dynamic_viscosity, medium.description());
        }

        auto const epsilon_T = medium
                                   .property(MaterialPropertyLib::PropertyType::
                                                 thermal_osmosis_permeability)
                                   .value<double>(variable_array, pos, t, dt);
        return epsilon_T * intrinsic_permeability / liquid_dynamic_viscosity;
    }

    if (medium.hasProperty(
            MaterialPropertyLib::PropertyType::thermal_osmosis_coefficient))
    {
        return MaterialPropertyLib::formEigenTensor<GlobalDim>(
            medium
                .property(MaterialPropertyLib::PropertyType::
                              thermal_osmosis_coefficient)
                .value(variable_array, pos, t, dt));
    }

    // Neither thermal_osmosis_permeability nor thermal_osmosis_coefficient is
    // defined.
    return Eigen::Matrix<double, GlobalDim, GlobalDim>::Zero();
}
}  // namespace ProcessLib
