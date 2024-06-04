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
    total_strain,
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
                             "total_strain",
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
    using Scalar = double;
    using KelvinVector = std::variant<std::monostate,
                                      Eigen::Matrix<double, 4, 1>,
                                      Eigen::Matrix<double, 6, 1>>;
    // Compare to GMatrixPolicy::GradientVectorType. The 1d case = Matrix<3, 1>
    // is not used so far.
    using DeformationGradient = std::variant<std::monostate,
                                             Eigen::Matrix<double, 5, 1>,
                                             Eigen::Matrix<double, 9, 1>>;

    using VariablePointerConst = std::
        variant<Scalar const*, KelvinVector const*, DeformationGradient const*>;

    VariablePointerConst address_of(Variable const v) const;

    using VariablePointer =
        std::variant<Scalar*, KelvinVector*, DeformationGradient*>;

    VariablePointer address_of(Variable const v);

    /// Read-only access.
    /// \note The returned value is a temporary.
    VariableType operator[](Variable const variable) const
    {
        return std::visit(
            []<typename T>(T* ptr) -> VariableType
            {
                auto identity = [](auto const& arg) -> VariableType
                { return arg; };
                if constexpr (std::is_same_v<Scalar const, T>)
                {
                    return *ptr;
                }
                else if constexpr (std::is_same_v<KelvinVector const, T>)
                {
                    return std::visit(identity, *ptr);
                }
                else if constexpr (std::is_same_v<DeformationGradient const, T>)
                {
                    return std::visit(identity, *ptr);
                }
                else
                {
                    static_assert(
                        !std::is_same_v<T, T>,
                        "Non-exhaustive visitor! The variable type (in the "
                        "std::is_same_v expression) must be one of the "
                        "VariableArray::{Scalar, KelvinVector, "
                        "DeformationGradient}.");
                }
            },
            address_of(variable));
    }

    double capillary_pressure = nan_;
    double concentration = nan_;
    DeformationGradient deformation_gradient;
    double density = nan_;
    double effective_pore_pressure = nan_;
    double enthalpy = nan_;
    double enthalpy_of_evaporation = nan_;
    double equivalent_plastic_strain = nan_;
    double grain_compressibility = nan_;
    double liquid_phase_pressure = nan_;
    double liquid_saturation = nan_;
    KelvinVector mechanical_strain;
    double molar_mass = nan_;
    double molar_mass_derivative = nan_;
    double molar_fraction = nan_;
    double gas_phase_pressure = nan_;
    double porosity = nan_;
    double solid_grain_pressure = nan_;
    KelvinVector stress;
    double temperature = nan_;
    KelvinVector total_strain;
    KelvinVector total_stress;
    double transport_porosity = nan_;
    double vapour_pressure = nan_;
    double volumetric_strain = nan_;

private:
    static constexpr auto nan_ = std::numeric_limits<double>::signaling_NaN();
};

static const VariableArray EmptyVariableArray{};

/// This function converts a string (e.g. a string from the configuration-tree)
/// into one of the entries of the Variable enumerator.
Variable convertStringToVariable(std::string const& string);
}  // namespace MaterialPropertyLib
