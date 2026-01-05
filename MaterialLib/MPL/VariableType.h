// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <BaseLib/Algorithm.h>
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
    fracture_aperture,
    grain_compressibility,
    ice_volume_fraction,
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
    volumetric_mechanical_strain,
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
                             "fracture_aperture",
                             "grain_compressibility",
                             "ice_volume_fraction",
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
                             "volumetric_mechanical_strain",
                             "volumetric_strain"}};

/// Data type for primary variables, designed to contain both scalar and vector
/// data.
using VariableType = std::variant<std::monostate,
                                  double,
                                  Eigen::Vector<double, 4>,
                                  Eigen::Vector<double, 5>,
                                  Eigen::Vector<double, 6>,
                                  Eigen::Vector<double, 9>>;

class VariableArray
{
public:
    using Scalar = double;
    using KelvinVector = std::variant<std::monostate,
                                      Eigen::Vector<double, 4>,
                                      Eigen::Vector<double, 6>>;
    // Compare to GMatrixPolicy::GradientVectorType. The 1d case = Vector<3>
    // is not used so far.
    using DeformationGradient = std::variant<std::monostate,
                                             Eigen::Vector<double, 5>,
                                             Eigen::Vector<double, 9>>;

    using VariablePointerConst = std::
        variant<Scalar const*, KelvinVector const*, DeformationGradient const*>;

    VariablePointerConst address_of(Variable const v) const;

    using VariablePointer =
        std::variant<Scalar*, KelvinVector*, DeformationGradient*>;

    template <typename Visitor>
    auto visitVariable(Visitor&& visitor, Variable const variable)
    {
        return std::visit(
            BaseLib::Overloaded{
                std::forward<Visitor>(visitor),
                []<typename T>(T*)
                {
                    static_assert(!std::is_same_v<T, T>,
                                  "Non-exhaustive visitor! The variable type "
                                  "must be one of the VariableArray::{Scalar, "
                                  "KelvinVector, DeformationGradient}.");
                }},
            address_of(variable));
    }

    template <typename Visitor>
    auto visitVariable(Visitor&& visitor, Variable const variable) const
    {
        return std::visit(
            BaseLib::Overloaded{
                std::forward<Visitor>(visitor),
                []<typename T>(T const*)
                {
                    static_assert(!std::is_same_v<T, T>,
                                  "Non-exhaustive visitor! The variable type "
                                  "must be one of the VariableArray::{Scalar, "
                                  "KelvinVector, DeformationGradient}.");
                }},
            address_of(variable));
    }

    /// Read-only access.
    /// \note The returned value is a temporary.
    VariableType operator[](Variable const variable) const
    {
        auto identity = [](auto const& arg) -> VariableType { return arg; };

        return visitVariable(
            BaseLib::Overloaded{
                [](Scalar const* ptr) -> VariableType { return *ptr; },
                [&identity](KelvinVector const* ptr) -> VariableType
                { return std::visit(identity, *ptr); },
                [&identity](DeformationGradient const* ptr) -> VariableType
                { return std::visit(identity, *ptr); }},
            variable);
    }

private:
    VariablePointer address_of(Variable const v);

public:
    double capillary_pressure = nan_;
    double concentration = nan_;
    DeformationGradient deformation_gradient;
    double density = nan_;
    double effective_pore_pressure = nan_;
    double enthalpy = nan_;
    double enthalpy_of_evaporation = nan_;
    double equivalent_plastic_strain = nan_;
    double fracture_aperture = nan_;
    double grain_compressibility = nan_;
    double ice_volume_fraction = nan_;
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
    double volumetric_mechanical_strain = nan_;
    double volumetric_strain = nan_;

    bool is2D() const;
    bool is3D() const;

private:
    static constexpr auto nan_ = std::numeric_limits<double>::signaling_NaN();
};

static const VariableArray EmptyVariableArray{};

/// This function converts a string (e.g. a string from the configuration-tree)
/// into one of the entries of the Variable enumerator.
Variable convertStringToVariable(std::string const& string);
}  // namespace MaterialPropertyLib
