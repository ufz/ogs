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

#include <Eigen/Dense>
#include <array>
#include <string>
#include <variant>

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
using SymmTensor = std::array<double, 6>;

/// Very simple 2d tensor data type for holding tensor components.
using Tensor2d = std::array<double, 4>;

/// Very simple tensor data type for holding
/// tensor components.
using Tensor = std::array<double, 9>;

/// Enum Variable is simply a list of all commonly used variables that are used
/// to determine the size of the VariableArray. If the variable of your choice
/// is missing, simply add it somewhere at the list, but above the last entry.
enum class Variable : int
{
    capillary_pressure,
    concentration,
    density,
    displacement,
    effective_pore_pressure,
    enthalpy_of_evaporation,
    equivalent_plastic_strain,
    grain_compressibility,
    liquid_phase_pressure,
    liquid_saturation,
    mechanical_strain,
    molar_mass,
    molar_fraction,
    phase_pressure,
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

/// Data type for primary variables, designed to contain both scalar and vector
/// data.
using VariableType =
    std::variant<std::monostate, double, Vector, Eigen::Matrix<double, 4, 1>,
                 Eigen::Matrix<double, 6, 1>>;

/// The VariableArray is a std::array of fixed size. Its size is determined by
/// the Variable enumerator list. Data type of that array is defined by the
/// VariableType definition.
using VariableArray =
    std::array<VariableType, static_cast<int>(Variable::number_of_variables)>;

Variable convertStringToVariable(std::string const& input);
}  // namespace MaterialPropertyLib
