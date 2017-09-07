/**
 * \author Norbert Grunwald
 * \date   Sep 7, 2017
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <boost/algorithm/string/predicate.hpp>
#include <boost/variant.hpp>
#include <array>
#include <cstddef>
#include <string>
#include "BaseLib/Error.h"

namespace MaterialPropertyLib
{
/// Very simple vector data type for holding
/// a pair of values.
using Pair = std::array<double, 2>;

/// Very simple vector data type for holding
/// vector components.
using Vector = std::array<double, 3>;

/// Simple symmetric tensor data type for holding
/// xx, yy, zz, xy, xz, yz tensor components.
using SymmTensor = std::array<double, 9>;

/// Very simple tensor data type for holding
/// tensor components.
using Tensor = std::array<double, 9>;

/// Variables is simply a list of all commonly used variables that are used to
/// determine the size of the VariableArray. If the variable of your choice is
/// missing, simply add it somewhere at the list, but above the last entry
/// 'numberOfVariables'
enum Variables : std::size_t
{
    phase_pressure,
    capillary_pressure,
    gas_density,
    liquid_density,
    temperature,
    liquid_saturation,
    u,
    numberOfVariables
};

/// Data type for primary variables, designed to contain both scalar and vector
/// data.
using VariableType = boost::variant<double, Vector>;

/// The VariableArray is a std::array of fixed size. Its size is determined by
/// the Variables enumerator list. Data type of that array is defined by the
/// VariableType definition.
using VariableArray = std::array<VariableType, numberOfVariables>;

/// This method returns a value of type double from the variables array
inline double getScalar(VariableType pv)
{
    return boost::get<double>(pv);
}

/// PropertyEnum is an enumerator list of all known properties of a substance.
/// This includes all properties on all scales (i.e. component, phase, and
/// medium scales). It is used as an indexer for the PropertyArray of the
/// materials. If a necessary property is not in the list, simply add a new one
/// in alphabetical order (of course, except for the last entry
/// 'number_of_property_enums'). Please note that any of these entries must also
/// appear in below convert functions.
enum PropertyEnum : std::size_t
{
    acentric_factor,
    binary_interaction_coefficient,
    biot_coefficient,
    brooks_corey_exponent,
    bulk_modulus,
    critical_density,
    critical_pressure,
    critical_temperature,
    compressibility,
    density,
    drhodT,
    effective_stress,
    entry_pressure,
    fredlund_parameters,
    heat_capacity,
    longitudinal_dispersivity,
    molar_mass,
    mole_fraction,
    molecular_diffusion,
    name,
    permeability,
    phase_velocity,
    porosity,
    reference_density,
    reference_temperature,
    reference_pressure,
    relative_permeability,
    residual_gas_saturation,
    residual_liquid_saturation,
    saturation,
    specific_heat_capacity,
    thermal_conductivity,
    thermal_expansivity,
    transveral_dispersivity,
    viscosity,
    number_of_property_enums
};

/// This function converts a string (e.g. a string from the configuration-tree)
/// into one of the entries of the PropertyEnum enumerator. To avoid confusion,
/// I suggest that the syntax of the properties in the input-files (i.e. the
/// strings) have to be identical to the syntax of the entries in the
/// enumerator.
inline PropertyEnum convertStringToProperty(std::string const& inString)
{
    if (boost::iequals(inString, "acentric_factor"))
    {
        return acentric_factor;
    }
    if (boost::iequals(inString, "binary_interaction_coefficient"))
    {
        return binary_interaction_coefficient;
    }
    if (boost::iequals(inString, "biot_coefficient"))
    {
        return biot_coefficient;
    }
    if (boost::iequals(inString, "brooks_corey_exponent"))
    {
        return brooks_corey_exponent;
    }
    if (boost::iequals(inString, "bulk_modulus"))
    {
        return bulk_modulus;
    }
    if (boost::iequals(inString, "critical_density"))
    {
        return critical_density;
    }
    if (boost::iequals(inString, "critical_pressure"))
    {
        return critical_pressure;
    }
    if (boost::iequals(inString, "critical_temperature"))
    {
        return critical_temperature;
    }
    if (boost::iequals(inString, "compressibility"))
    {
        return compressibility;
    }
    if (boost::iequals(inString, "density"))
    {
        return density;
    }
    if (boost::iequals(inString, "drhodT"))
    {
        return drhodT;
    }
    if (boost::iequals(inString, "effective_stress"))
    {
        return effective_stress;
    }
    if (boost::iequals(inString, "entry_pressure"))
    {
        return entry_pressure;
    }
    if (boost::iequals(inString, "fredlund_parameters"))
    {
        return fredlund_parameters;
    }
    if (boost::iequals(inString, "heat_capacity"))
    {
        return heat_capacity;
    }
    if (boost::iequals(inString, "longitudinal_dispersivity"))
    {
        return longitudinal_dispersivity;
    }
    if (boost::iequals(inString, "molar_mass"))
    {
        return molar_mass;
    }
    if (boost::iequals(inString, "mole_fraction"))
    {
        return mole_fraction;
    }
    if (boost::iequals(inString, "molecular_diffusion"))
    {
        return molecular_diffusion;
    }
    if (boost::iequals(inString, "name"))
    {
        return name;
    }
    if (boost::iequals(inString, "permeability"))
    {
        return permeability;
    }
    if (boost::iequals(inString, "porosity"))
    {
        return porosity;
    }
    if (boost::iequals(inString, "phase_velocity"))
    {
        return phase_velocity;
    }
    if (boost::iequals(inString, "reference_density"))
    {
        return reference_density;
    }
    if (boost::iequals(inString, "reference_temperature"))
    {
        return reference_temperature;
    }
    if (boost::iequals(inString, "reference_pressure"))
    {
        return reference_pressure;
    }
    if (boost::iequals(inString, "relative_permeability"))
    {
        return relative_permeability;
    }
    if (boost::iequals(inString, "residual_gas_saturation"))
    {
        return residual_gas_saturation;
    }
    if (boost::iequals(inString, "residual_liquid_saturation"))
    {
        return residual_liquid_saturation;
    }
    if (boost::iequals(inString, "saturation"))
    {
        return saturation;
    }
    if (boost::iequals(inString, "specific_heat_capacity"))
    {
        return specific_heat_capacity;
    }
    if (boost::iequals(inString, "thermal_conductivity"))
    {
        return thermal_conductivity;
    }
    if (boost::iequals(inString, "thermal_expansivity"))
    {
        return thermal_expansivity;
    }
    if (boost::iequals(inString, "transveral_dispersivity"))
    {
        return transveral_dispersivity;
    }
    if (boost::iequals(inString, "viscosity"))
    {
        return viscosity;
    }

    OGS_FATAL(
        "The property name \"%s\" does not correspond to any known "
        "property",
        inString.c_str());

    return number_of_property_enums;  // to avoid the 'no return' warning
}

const static std::vector<std::string> convertEnumToString{
    {"acentric_factor"},
    {"binary_interaction_coefficient"},
    {"biot_coefficient"},
    {"brooks_corey_exponent"},
    {"bulk_modulus"},
    {"critical_density"},
    {"critical_pressure"},
    {"critical_temperature"},
    {"compressibility"},
    {"density"},
    {"drhodT"},
    {"effective_stress"},
    {"entry_pressure"},
    {"fredlund_parameters"},
    {"heat_capacity"},
    {"longitudinal_dispersivity"},
    {"molar_mass"},
    {"mole_fraction"},
    {"molecular_diffusion"},
    {"name"},
    {"permeability"},
    {"phase_velocity"},
    {"porosity"},
    {"reference_density"},
    {"reference_temperature"},
    {"reference_pressure"},
    {"relative_permeability"},
    {"residual_gas_saturation"},
    {"residual_liquid_saturation"},
    {"saturation"},
    {"specific_heat_capacity"},
    {"thermal_conductivity"},
    {"thermal_expansivity"},
    {"transveral_dispersivity"},
    {"viscosity"}};

}  // namespace MaterialPropertyLib
