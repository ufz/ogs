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

#include "BaseLib/Error.h"
#include <boost/algorithm/string/predicate.hpp>
#include <cstddef>
#include <string>

namespace MaterialPropertyLib
{

enum PropertyEnum : std::size_t
{
    molar_mass,
    critical_temperature,
    density,
    viscosity,
    thermal_conductivity,
    effective_stress,
    permeability,
    phase_velocity,
    saturation,
    relative_permeability,
    number_of_property_enums
};

inline PropertyEnum convertStringToProperty (std::string const& inString) {
    if (boost::iequals(inString, "molar_mass")) return molar_mass;
    if (boost::iequals(inString, "critical_temperature")) return critical_temperature;
    if (boost::iequals(inString, "density")) return density;
    if (boost::iequals(inString, "viscosity")) return viscosity;
    if (boost::iequals(inString, "thermal_conductivity")) return thermal_conductivity;
    if (boost::iequals(inString, "effective_stress")) return effective_stress;
    if (boost::iequals(inString, "permeability")) return permeability;
    if (boost::iequals(inString, "phase_velocity")) return phase_velocity;
    if (boost::iequals(inString, "saturation")) return saturation;
    if (boost::iequals(inString, "relative_permeability")) return relative_permeability;

    OGS_FATAL("The property name \"%s\" does not correspond to any known "
            "property", inString.c_str());
    return number_of_property_enums;  // to avoid the 'no return' warning
}



}   //MaterialPropertyLib

#endif /* MATERIALLIB_MPL_MPENUMS_H_ */
