// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "BaseLib/Error.h"
#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Phase.h"
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
/// \f$k_T = \epsilon_T k / \mu\f$ with the thermo-osmotic permeability
/// \f$\epsilon_T\f$ (\f$[\epsilon_T] = Pa/K\f$), which converts a temperature
/// gradient into the equivalent pore pressure gradient driving the
/// thermo-osmotic flux, the intrinsic permeability
/// \f$k\f$ (\f$[k] = m^2\f$), and the liquid's dynamic viscosity
/// \f$\mu\f$ (\f$[\mu] = Pa \cdot s\f$).
/// Defining both properties at the same time is an error. If neither is
/// defined, the zero tensor is returned (no thermo-osmosis).
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
/// \note The product below is formed as \f$\epsilon_T k\f$. That order is
/// only exercised for diagonal \f$\epsilon_T\f$ and \f$k\f$, for which the
/// two tensors commute and the order does not matter. For non-diagonal
/// tensors the order matters, is not covered by any test, and may need to be
/// swapped if it turns out not to be the physically intended composition.
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
    auto const solid_phase =
        getOptionalPhase(medium, MaterialPropertyLib::PhaseName::Solid);
    if (solid_phase &&
        solid_phase->hasProperty(
            MaterialPropertyLib::PropertyType::thermal_osmosis_coefficient))
    {
        OGS_FATAL(
            "thermal_osmosis_coefficient is defined on the solid phase of "
            "{:s}, but is now read from the medium. Move the property from "
            "the solid phase to the medium, or use "
            "scripts/dev/move_thermal_osmosis.py to migrate the project "
            "file.",
            medium.description());
    }

    bool const has_thermal_osmosis_permeability = medium.hasProperty(
        MaterialPropertyLib::PropertyType::thermal_osmosis_permeability);

    bool const has_thermal_osmosis_coefficient = medium.hasProperty(
        MaterialPropertyLib::PropertyType::thermal_osmosis_coefficient);

    if (has_thermal_osmosis_permeability && has_thermal_osmosis_coefficient)
    {
        OGS_FATAL(
            "Thermo-osmosis permeability and coefficient cannot be defined at "
            "the same time in {:s}.",
            medium.description());
    }

    if (has_thermal_osmosis_permeability)
    {
        if (liquid_dynamic_viscosity <= 0.)
        {
            OGS_FATAL(
                "Liquid dynamic viscosity must be > 0 when using "
                "thermal_osmosis_permeability, but is {:g} in {:s}.",
                liquid_dynamic_viscosity, medium.description());
        }

        auto const epsilon_T = MaterialPropertyLib::formEigenTensor<GlobalDim>(
            medium
                .property(MaterialPropertyLib::PropertyType::
                              thermal_osmosis_permeability)
                .value(variable_array, pos, t, dt));
        return epsilon_T * intrinsic_permeability / liquid_dynamic_viscosity;
    }

    if (has_thermal_osmosis_coefficient)
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
