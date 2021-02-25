/**
 * \file
 * \author Norbert Grunwald
 * \date   Sep 7, 2017
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <array>
#include <boost/algorithm/string/predicate.hpp>
#include <memory>
#include <string>

#include "BaseLib/Error.h"

namespace MaterialPropertyLib
{
class Property;
}

namespace MaterialPropertyLib
{
/// PropertyType is an enumerator list of all known properties of a substance.
/// This includes all properties on all scales (i.e. component, phase, and
/// medium scales). It is used as an index for the PropertyArray of the
/// materials. If a necessary property is not in the list, simply add a new one
/// in alphabetical order (of course, except for the last entry). Please note
/// that any of these entries must also appear in below convert functions.
enum PropertyType : int
{
    acentric_factor,
    binary_interaction_coefficient,
    biot_coefficient,
    bishops_effective_stress,
    brooks_corey_exponent,
    bulk_modulus,
    capillary_pressure,
    critical_density,
    critical_pressure,
    critical_temperature,
    compressibility,
    /// used to specify decay rate of a substance.
    decay_rate,
    density,
    diffusion,
    drhodT,
    effective_stress,
    entry_pressure,
    evaporation_enthalpy,
    fredlund_parameters,
    heat_capacity,
    /// used to compute the hydrodynamic dispersion tensor.
    longitudinal_dispersivity,
    molality,
    molar_mass,
    molar_volume,
    mole_fraction,
    /// ion diffusivity in free water.
    molecular_diffusion,
    name,
    permeability,
    phase_velocity,
    /// ion diffusivity in the porous medium with account of the effect of
    /// tortuosity and connectivity.
    pore_diffusion,
    porosity,
    reference_density,
    reference_temperature,
    reference_pressure,
    relative_permeability,
    relative_permeability_nonwetting_phase,
    residual_gas_saturation,
    residual_liquid_saturation,
    /// specify retardation factor used in component transport process.
    retardation_factor,
    saturation,
    /// capillary pressure saturation relationship for microstructure.
    saturation_micro,
    specific_heat_capacity,
    storage,
    swelling_stress_rate,
    thermal_conductivity,
    thermal_expansivity,
    thermal_longitudinal_dispersivity,
    thermal_osmosis_coefficient,
    thermal_transversal_dispersivity,
    transport_porosity,
    /// used to compute the hydrodynamic dispersion tensor.
    transversal_dispersivity,
    vapor_pressure,
    viscosity,
    volume_fraction,
    number_of_properties
};

/// This function converts a string (e.g. a string from the configuration-tree)
/// into one of the entries of the PropertyType enumerator. To avoid confusion,
/// I suggest that the syntax of the properties in the input-files (i.e. the
/// strings) have to be identical to the syntax of the entries in the
/// enumerator.
inline PropertyType convertStringToProperty(std::string const& inString)
{
    if (boost::iequals(inString, "acentric_factor"))
    {
        return PropertyType::acentric_factor;
    }
    if (boost::iequals(inString, "binary_interaction_coefficient"))
    {
        return PropertyType::binary_interaction_coefficient;
    }
    if (boost::iequals(inString, "biot_coefficient"))
    {
        return PropertyType::biot_coefficient;
    }
    if (boost::iequals(inString, "bishops_effective_stress"))
    {
        return PropertyType::bishops_effective_stress;
    }
    if (boost::iequals(inString, "brooks_corey_exponent"))
    {
        return PropertyType::brooks_corey_exponent;
    }
    if (boost::iequals(inString, "bulk_modulus"))
    {
        return PropertyType::bulk_modulus;
    }
    if (boost::iequals(inString, "capillary_pressure"))
    {
        return PropertyType::capillary_pressure;
    }
    if (boost::iequals(inString, "critical_density"))
    {
        return PropertyType::critical_density;
    }
    if (boost::iequals(inString, "critical_pressure"))
    {
        return PropertyType::critical_pressure;
    }
    if (boost::iequals(inString, "critical_temperature"))
    {
        return PropertyType::critical_temperature;
    }
    if (boost::iequals(inString, "compressibility"))
    {
        return PropertyType::compressibility;
    }
    if (boost::iequals(inString, "decay_rate"))
    {
        return PropertyType::decay_rate;
    }
    if (boost::iequals(inString, "density"))
    {
        return PropertyType::density;
    }
    if (boost::iequals(inString, "diffusion"))
    {
        return PropertyType::diffusion;
    }
    if (boost::iequals(inString, "drhodT"))
    {
        return PropertyType::drhodT;
    }
    if (boost::iequals(inString, "effective_stress"))
    {
        return PropertyType::effective_stress;
    }
    if (boost::iequals(inString, "entry_pressure"))
    {
        return PropertyType::entry_pressure;
    }
    if (boost::iequals(inString, "evaporation_enthalpy"))
    {
        return PropertyType::evaporation_enthalpy;
    }
    if (boost::iequals(inString, "fredlund_parameters"))
    {
        return PropertyType::fredlund_parameters;
    }
    if (boost::iequals(inString, "heat_capacity"))
    {
        return PropertyType::heat_capacity;
    }
    if (boost::iequals(inString, "longitudinal_dispersivity"))
    {
        return PropertyType::longitudinal_dispersivity;
    }
    if (boost::iequals(inString, "molality"))
    {
        return PropertyType::molality;
    }
    if (boost::iequals(inString, "molar_mass"))
    {
        return PropertyType::molar_mass;
    }
    if (boost::iequals(inString, "molar_volume"))
    {
        return PropertyType::molar_volume;
    }
    if (boost::iequals(inString, "mole_fraction"))
    {
        return PropertyType::mole_fraction;
    }
    if (boost::iequals(inString, "molecular_diffusion"))
    {
        return PropertyType::molecular_diffusion;
    }
    if (boost::iequals(inString, "name"))
    {
        return PropertyType::name;
    }
    if (boost::iequals(inString, "permeability"))
    {
        return PropertyType::permeability;
    }
    if (boost::iequals(inString, "pore_diffusion"))
    {
        return PropertyType::pore_diffusion;
    }
    if (boost::iequals(inString, "porosity"))
    {
        return PropertyType::porosity;
    }
    if (boost::iequals(inString, "phase_velocity"))
    {
        return PropertyType::phase_velocity;
    }
    if (boost::iequals(inString, "reference_density"))
    {
        return PropertyType::reference_density;
    }
    if (boost::iequals(inString, "reference_temperature"))
    {
        return PropertyType::reference_temperature;
    }
    if (boost::iequals(inString, "reference_pressure"))
    {
        return PropertyType::reference_pressure;
    }
    if (boost::iequals(inString, "relative_permeability"))
    {
        return PropertyType::relative_permeability;
    }
    if (boost::iequals(inString, "relative_permeability_nonwetting_phase"))
    {
        return PropertyType::relative_permeability_nonwetting_phase;
    }
    if (boost::iequals(inString, "residual_gas_saturation"))
    {
        return PropertyType::residual_gas_saturation;
    }
    if (boost::iequals(inString, "residual_liquid_saturation"))
    {
        return PropertyType::residual_liquid_saturation;
    }
    if (boost::iequals(inString, "retardation_factor"))
    {
        return PropertyType::retardation_factor;
    }
    if (boost::iequals(inString, "saturation"))
    {
        return PropertyType::saturation;
    }
    if (boost::iequals(inString, "saturation_micro"))
    {
        return PropertyType::saturation_micro;
    }
    if (boost::iequals(inString, "specific_heat_capacity"))
    {
        return PropertyType::specific_heat_capacity;
    }
    if (boost::iequals(inString, "storage"))
    {
        return PropertyType::storage;
    }
    if (boost::iequals(inString, "swelling_stress_rate"))
    {
        return PropertyType::swelling_stress_rate;
    }
    if (boost::iequals(inString, "thermal_conductivity"))
    {
        return PropertyType::thermal_conductivity;
    }
    if (boost::iequals(inString, "thermal_expansivity"))
    {
        return PropertyType::thermal_expansivity;
    }
    if (boost::iequals(inString, "thermal_longitudinal_dispersivity"))
    {
        return PropertyType::thermal_longitudinal_dispersivity;
    }
    if (boost::iequals(inString, "thermal_osmosis_coefficient"))
    {
        return PropertyType::thermal_osmosis_coefficient;
    }
    if (boost::iequals(inString, "thermal_transversal_dispersivity"))
    {
        return PropertyType::thermal_transversal_dispersivity;
    }
    if (boost::iequals(inString, "transport_porosity"))
    {
        return PropertyType::transport_porosity;
    }
    if (boost::iequals(inString, "transversal_dispersivity"))
    {
        return PropertyType::transversal_dispersivity;
    }
    if (boost::iequals(inString, "vapor_pressure"))
    {
        return PropertyType::vapor_pressure;
    }
    if (boost::iequals(inString, "viscosity"))
    {
        return PropertyType::viscosity;
    }
    if (boost::iequals(inString, "volume_fraction"))
    {
        return PropertyType::volume_fraction;
    }

    OGS_FATAL(
        "The property name '{:s}' does not correspond to any known property",
        inString);

    return PropertyType::number_of_properties;  // to avoid the 'no return'
                                                // warning
}

static const std::array<std::string, PropertyType::number_of_properties>
    property_enum_to_string{{"acentric_factor",
                             "binary_interaction_coefficient",
                             "biot_coefficient",
                             "bishops_effective_stress",
                             "brooks_corey_exponent",
                             "bulk_modulus",
                             "capillary_pressure",
                             "critical_density",
                             "critical_pressure",
                             "critical_temperature",
                             "compressibility",
                             "decay_rate",
                             "density",
                             "diffusion",
                             "drhodT",
                             "effective_stress",
                             "entry_pressure",
                             "evaporation_enthalpy",
                             "fredlund_parameters",
                             "heat_capacity",
                             "longitudinal_dispersivity",
                             "molality",
                             "molar_mass",
                             "molar_volume",
                             "mole_fraction",
                             "molecular_diffusion",
                             "name",
                             "permeability",
                             "phase_velocity",
                             "pore_diffusion",
                             "porosity",
                             "reference_density",
                             "reference_temperature",
                             "reference_pressure",
                             "relative_permeability",
                             "relative_permeability_nonwetting_phase",
                             "residual_gas_saturation",
                             "residual_liquid_saturation",
                             "retardation_factor",
                             "saturation",
                             "saturation_micro",
                             "specific_heat_capacity",
                             "storage",
                             "swelling_stress_rate",
                             "thermal_conductivity",
                             "thermal_expansivity",
                             "thermal_longitudinal_dispersivity",
                             "thermal_osmosis_coefficient",
                             "thermal_transversal_dispersivity",
                             "transport_porosity",
                             "transversal_dispersivity",
                             "vapor_pressure",
                             "viscosity",
                             "volume_fraction"}};

/// This data type is based on a std::array. It can hold pointers to objects of
/// class Property or its inheritors. The size of this array is determined by
/// the number of entries of the PropertyType enumerator.
using PropertyArray =
    std::array<std::unique_ptr<Property>, PropertyType::number_of_properties>;

}  // namespace MaterialPropertyLib
