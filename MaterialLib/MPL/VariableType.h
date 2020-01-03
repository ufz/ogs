/**
 * \file
 * \author Norbert Grunwald
 * \date   Sep 7, 2017
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <array>
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
    concentration,
    phase_pressure,
    capillary_pressure,
    density,
    temperature,
    liquid_saturation,
    displacement,
    number_of_variables
};

/// Data type for primary variables, designed to contain both scalar and vector
/// data.
using VariableType = std::variant<double, Vector>;

/// The VariableArray is a std::array of fixed size. Its size is determined by
/// the Variable enumerator list. Data type of that array is defined by the
/// VariableType definition.
using VariableArray =
    std::array<VariableType, static_cast<int>(Variable::number_of_variables)>;

/// This method returns a value of type double from the variables array
inline double getScalar(VariableType pv)
{
    return std::get<double>(pv);
}

Variable convertStringToVariable(std::string const& input);
}  // namespace MaterialPropertyLib
