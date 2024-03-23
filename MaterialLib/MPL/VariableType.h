/**
 * \file
 * \author Norbert Grunwald
 * \date   Sep 7, 2017
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <BaseLib/Error.h>

#include <Eigen/Core>
#include <array>
#include <string>
#include <variant>

namespace MaterialPropertyLib
{

/// Enum Variable is simply a list of all commonly used variables.
/// If the variable of your choice is missing, simply add it somewhere at the
/// list, but above the last entry.
enum class Variable : int
{
    capillary_pressure,
    concentration,
    deformation_gradient,
    density,
    effective_pore_pressure,
    enthalpy,
    enthalpy_of_evaporation,
    equivalent_plastic_strain,
    grain_compressibility,
    liquid_phase_pressure,
    liquid_saturation,
    mechanical_strain,
    molar_mass,
    molar_mass_derivative,
    molar_fraction,
    gas_phase_pressure,
    porosity,
    solid_grain_pressure,
    stress,
    temperature,
    total_stress,
    transport_porosity,
    vapour_pressure,
    volumetric_strain,
    number_of_variables
};

static const std::array<std::string,
                        static_cast<int>(Variable::number_of_variables)>
    variable_enum_to_string{{"capillary_pressure",
                             "concentration",
                             "deformation_gradient",
                             "density",
                             "effective_pore_pressure",
                             "enthalpy",
                             "enthalpy_of_evaporation",
                             "equivalent_plastic_strain",
                             "grain_compressibility",
                             "liquid_phase_pressure",
                             "liquid_saturation",
                             "mechanical_strain",
                             "molar_mass",
                             "molar_mass_derivative",
                             "molar_fraction",
                             "gas_phase_pressure",
                             "porosity",
                             "solid_grain_pressure",
                             "stress",
                             "temperature",
                             "total_stress",
                             "transport_porosity",
                             "vapour_pressure",
                             "volumetric_strain"}};

/// Data type for primary variables, designed to contain both scalar and vector
/// data.
using VariableType = std::variant<std::monostate,
                                  double,
                                  Eigen::Matrix<double, 4, 1>,
                                  Eigen::Matrix<double, 5, 1>,
                                  Eigen::Matrix<double, 6, 1>,
                                  Eigen::Matrix<double, 9, 1>>;

class VariableArray
{
public:
    /// Read-only access.
    /// \note The returned value is a temporary.
    VariableType operator[](Variable const variable) const
    {
        auto identity = [](auto&& arg) -> VariableType { return arg; };
        switch (variable)
        {
            case Variable::capillary_pressure:
                return capillary_pressure;
            case Variable::concentration:
                return concentration;
            case Variable::deformation_gradient:
                return std::visit(identity, deformation_gradient);
            case Variable::density:
                return density;
            case Variable::effective_pore_pressure:
                return effective_pore_pressure;
            case Variable::enthalpy:
                return enthalpy;
            case Variable::enthalpy_of_evaporation:
                return enthalpy_of_evaporation;
            case Variable::equivalent_plastic_strain:
                return equivalent_plastic_strain;
            case Variable::grain_compressibility:
                return grain_compressibility;
            case Variable::liquid_phase_pressure:
                return liquid_phase_pressure;
            case Variable::liquid_saturation:
                return liquid_saturation;
            case Variable::mechanical_strain:
                return std::visit(identity, mechanical_strain);
            case Variable::molar_mass:
                return molar_mass;
            case Variable::molar_mass_derivative:
                return molar_mass_derivative;
            case Variable::molar_fraction:
                return molar_fraction;
            case Variable::gas_phase_pressure:
                return gas_phase_pressure;
            case Variable::porosity:
                return porosity;
            case Variable::solid_grain_pressure:
                return solid_grain_pressure;
            case Variable::stress:
                return std::visit(identity, stress);
            case Variable::temperature:
                return temperature;
            case Variable::total_stress:
                return std::visit(identity, total_stress);
            case Variable::transport_porosity:
                return transport_porosity;
            case Variable::vapour_pressure:
                return vapour_pressure;
            case Variable::volumetric_strain:
                return volumetric_strain;
            default:
                OGS_FATAL(
                    "No conversion to VariableType is provided for variable "
                    "{:d}",
                    static_cast<int>(variable));
        };
    }

    double capillary_pressure = nan_;
    double concentration = nan_;
    // Compare to GMatrixPolicy::GradientVectorType. The 1d case = Matrix<3, 1>
    // is not used so far.
    std::variant<std::monostate,
                 Eigen::Matrix<double, 5, 1>,
                 Eigen::Matrix<double, 9, 1>>
        deformation_gradient;
    double density = nan_;
    double effective_pore_pressure = nan_;
    double enthalpy = nan_;
    double enthalpy_of_evaporation = nan_;
    double equivalent_plastic_strain = nan_;
    double grain_compressibility = nan_;
    double liquid_phase_pressure = nan_;
    double liquid_saturation = nan_;
    std::variant<std::monostate,
                 Eigen::Matrix<double, 4, 1>,
                 Eigen::Matrix<double, 6, 1>>
        mechanical_strain;
    double molar_mass = nan_;
    double molar_mass_derivative = nan_;
    double molar_fraction = nan_;
    double gas_phase_pressure = nan_;
    double porosity = nan_;
    double solid_grain_pressure = nan_;
    std::variant<std::monostate,
                 Eigen::Matrix<double, 4, 1>,
                 Eigen::Matrix<double, 6, 1>>
        stress;
    double temperature = nan_;
    std::variant<std::monostate,
                 Eigen::Matrix<double, 4, 1>,
                 Eigen::Matrix<double, 6, 1>>
        total_strain;
    std::variant<std::monostate,
                 Eigen::Matrix<double, 4, 1>,
                 Eigen::Matrix<double, 6, 1>>
        total_stress;
    double transport_porosity = nan_;
    double vapour_pressure = nan_;
    double volumetric_strain = nan_;

private:
    static constexpr auto nan_ = std::numeric_limits<double>::signaling_NaN();
};

static const VariableArray EmptyVariableArray{};

/// This function converts a string (e.g. a string from the configuration-tree)
/// into one of the entries of the VariableType enumerator.
Variable convertStringToVariable(std::string const& string);
}  // namespace MaterialPropertyLib
