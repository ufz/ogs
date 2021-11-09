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
#include <memory>
#include <string>

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
    concentration,
    decay_rate,
    density,
    diffusion,
    drhodT,
    effective_stress,
    entry_pressure,
    evaporation_enthalpy,
    fredlund_parameters,
    heat_capacity,
    latent_heat,
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
    poissons_ratio,
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
    specific_latent_heat,
    storage,
    storage_contribution,
    swelling_stress_rate,
    thermal_conductivity,
    /// Thermal diffusion enhancement factor for water vapor flow
    thermal_diffusion_enhancement_factor,
    /// The thermal expansivity corresponds to the linear thermal expansion
    /// coefficient for a solid and to the volumetric thermal expansion
    /// coefficient for a fluid
    thermal_expansivity,
    thermal_expansivity_contribution,
    thermal_longitudinal_dispersivity,
    thermal_osmosis_coefficient,
    thermal_transversal_dispersivity,
    transport_porosity,
    /// used to compute the hydrodynamic dispersion tensor.
    transversal_dispersivity,
    vapour_pressure,
    vapour_density,
    vapour_diffusion,
    viscosity,
    volume_fraction,
    youngs_modulus,
    number_of_properties
};

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
                             "concentration",
                             "decay_rate",
                             "density",
                             "diffusion",
                             "drhodT",
                             "effective_stress",
                             "entry_pressure",
                             "evaporation_enthalpy",
                             "fredlund_parameters",
                             "heat_capacity",
                             "latent_heat",
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
                             "poissons_ratio",
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
                             "specific_latent_heat",
                             "storage",
                             "storage_contribution",
                             "swelling_stress_rate",
                             "thermal_conductivity",
                             "thermal_diffusion_enhancement_factor",
                             "thermal_expansivity",
                             "thermal_expansivity_contribution",
                             "thermal_longitudinal_dispersivity",
                             "thermal_osmosis_coefficient",
                             "thermal_transversal_dispersivity",
                             "transport_porosity",
                             "transversal_dispersivity",
                             "vapour_pressure",
                             "vapour_density",
                             "vapour_diffusion",
                             "viscosity",
                             "volume_fraction",
                             "youngs_modulus"}};

/// This function converts a string (e.g. a string from the configuration-tree)
/// into one of the entries of the PropertyType enumerator.
PropertyType convertStringToProperty(std::string const& string);

/// This data type is based on a std::array. It can hold pointers to objects of
/// class Property or its inheritors. The size of this array is determined by
/// the number of entries of the PropertyType enumerator.
using PropertyArray =
    std::array<std::unique_ptr<Property>, PropertyType::number_of_properties>;

}  // namespace MaterialPropertyLib
