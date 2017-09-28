/**
 * \author Norbert Grunwald
 * \date   Sep 7, 2017
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MATERIALLIB_MPL_MPENUMS_H_
#define MATERIALLIB_MPL_MPENUMS_H_

#include <boost/algorithm/string/predicate.hpp>
#include <boost/variant.hpp>
#include <cstddef>
#include <string>
#include "BaseLib/Error.h"

namespace MaterialPropertyLib
{
/// Very simple vector data type for holding
/// vector components.
using Vector = std::array<double, 3>;
/// Very simple tensor data type for holding
/// tensor components.
using Tensor = std::array<double, 9>;

/**
 * PrimaryVariables is simply a list of all commonly used primary variables.
 * It is used to determine the size of the VariableArray. If the primary
 * variable of your choice is missing, simply add it somewhere at the list,
 * but above the last entry 'numberOfPrimaryVariables'
*/
enum PrimaryVariables : std::size_t
{
    p_GR,
    p_LR,
    p_cap,
    T_S,
    T_L,
    T_G,
    T,
    u,
    numberOfPrimaryVariables
};

/**
 * Data type for primary variables, designed to contain both
 * scalar and vector data.
 */
using PrimaryVariableType = boost::variant<double, Vector>;

/**
 * The VariableArray is a std::array of fixed size. Its size is determined
 * by the PrimaryVariables enumerator list. Data type of that array is
 * defined by the PrimaryVariableType definition.
*/
using VariableArray = std::array<PrimaryVariableType, numberOfPrimaryVariables>;

/// This method returns a value of type double from the
/// primar variables array
inline double getScalar(PrimaryVariableType pv)
{
    return boost::get<double>(pv);
}

/**
 * PropertyEnum is an enumerator list of all known properties of a substance.
 * This includes all properties on all scales (i.e. component, phase, and
 * medium scales). It is used as an indexer for the PropertyArray of the
 * materials. If a necessary property is not in the list, simply add a new
 * one in alphabetical order (of course, except for the last entry
 * 'number_of_property_enums'). Please note that any of these entries must
 * also appear in below convert functions.
*/
enum PropertyEnum : std::size_t
{
    acentric_factor,
    binary_interaction_coefficient,
    critical_density,
    critical_pressure,
    critical_temperature,
    density,
    drhodT,
    effective_stress,
    heat_capacity,
    molar_mass,
    mole_fraction,
    name,
    permeability,
    phase_velocity,
    reference_density,
    reference_temperature,
    relative_permeability,
    saturation,
    thermal_conductivity,
    viscosity,
    number_of_property_enums
};

/**
 * This function converts a string (e.g. a string from the configuration-tree)
 * into one of the entries of the PropertyEnum enumerator. To avoid confusion,
 * I suggest that the syntax of the properties in the input-files (i.e. the
 * strings) have to be identical to the syntax of the entries in the
 * enumerator.
*/
inline PropertyEnum convertStringToProperty(std::string const& inString)
{
    if (boost::iequals(inString, "acentric_factor"))
        return acentric_factor;
    if (boost::iequals(inString, "binary_interaction_coefficient"))
        return binary_interaction_coefficient;
    if (boost::iequals(inString, "critical_density"))
        return critical_density;
    if (boost::iequals(inString, "critical_pressure"))
        return critical_pressure;
    if (boost::iequals(inString, "critical_temperature"))
        return critical_temperature;
    if (boost::iequals(inString, "density"))
        return density;
    if (boost::iequals(inString, "drhodT"))
        return drhodT;
    if (boost::iequals(inString, "effective_stress"))
        return effective_stress;
    if (boost::iequals(inString, "heat_capacity"))
        return heat_capacity;
    if (boost::iequals(inString, "molar_mass"))
        return molar_mass;
    if (boost::iequals(inString, "mole_fraction"))
        return mole_fraction;
    if (boost::iequals(inString, "name"))
        return name;
    if (boost::iequals(inString, "permeability"))
        return permeability;
    if (boost::iequals(inString, "phase_velocity"))
        return phase_velocity;
    if (boost::iequals(inString, "reference_density"))
        return reference_density;
    if (boost::iequals(inString, "reference_temperature"))
        return reference_temperature;
    if (boost::iequals(inString, "relative_permeability"))
        return relative_permeability;
    if (boost::iequals(inString, "saturation"))
        return saturation;
    if (boost::iequals(inString, "thermal_conductivity"))
        return thermal_conductivity;
    if (boost::iequals(inString, "viscosity"))
        return viscosity;

    OGS_FATAL(
        "The property name \"%s\" does not correspond to any known "
        "property",
        inString.c_str());

    return number_of_property_enums;  // to avoid the 'no return' warning
}

const static std::vector<std::string> convertEnumToString{
    {"acentric_factor"},
    {"binary_interaction_coefficient"},
    {"critical_density"},
    {"critical_pressure"},
    {"critical_temperature"},
    {"density"},
    {"drhodT"},
    {"effective_stress"},
    {"heat_capacity"},
    {"molar_mass"},
    {"mole_fraction"},
    {"name"},
    {"permeability"},
    {"phase_velocity"},
    {"reference_density"},
    {"reference_temperature"},
    {"relative_permeability"},
    {"saturation"},
    {"thermal_conductivity"},
    {"viscosity"}};

}  // MaterialPropertyLib

#endif /* MATERIALLIB_MPL_MPENUMS_H_ */
