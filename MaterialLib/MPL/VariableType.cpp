/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "VariableType.h"

#include <boost/algorithm/string/predicate.hpp>

#include "BaseLib/Error.h"

/// Returns true if the argument is uninitialized or has Rows number of rows.
template <int Rows>
static bool maybeHasSize(const auto& arg)
{
    return std::holds_alternative<std::monostate>(arg) ||
           std::holds_alternative<Eigen::Matrix<double, Rows, 1>>(arg);
}

namespace MaterialPropertyLib
{
Variable convertStringToVariable(std::string const& string)
{
    for (int i = 0; i < static_cast<int>(Variable::number_of_variables); ++i)
    {
        if (boost::iequals(string, variable_enum_to_string[i]))
        {
            return static_cast<Variable>(i);
        }
    }

    OGS_FATAL(
        "The variable name '{:s}' does not correspond to any known variable",
        string);
}

VariableArray::VariablePointerConst VariableArray::address_of(
    Variable const v) const
{
    switch (v)
    {
        case Variable::capillary_pressure:
            return &capillary_pressure;
        case Variable::concentration:
            return &concentration;
        case Variable::deformation_gradient:
            return &deformation_gradient;
        case Variable::density:
            return &density;
        case Variable::effective_pore_pressure:
            return &effective_pore_pressure;
        case Variable::enthalpy:
            return &enthalpy;
        case Variable::enthalpy_of_evaporation:
            return &enthalpy_of_evaporation;
        case Variable::equivalent_plastic_strain:
            return &equivalent_plastic_strain;
        case Variable::fracture_aperture:
            return &fracture_aperture;
        case Variable::grain_compressibility:
            return &grain_compressibility;
        case Variable::ice_volume_fraction:
            return &ice_volume_fraction;
        case Variable::liquid_phase_pressure:
            return &liquid_phase_pressure;
        case Variable::liquid_saturation:
            return &liquid_saturation;
        case Variable::mechanical_strain:
            return &mechanical_strain;
        case Variable::molar_mass:
            return &molar_mass;
        case Variable::molar_mass_derivative:
            return &molar_mass_derivative;
        case Variable::molar_fraction:
            return &molar_fraction;
        case Variable::gas_phase_pressure:
            return &gas_phase_pressure;
        case Variable::porosity:
            return &porosity;
        case Variable::solid_grain_pressure:
            return &solid_grain_pressure;
        case Variable::stress:
            return &stress;
        case Variable::temperature:
            return &temperature;
        case Variable::total_strain:
            return &total_strain;
        case Variable::total_stress:
            return &total_stress;
        case Variable::transport_porosity:
            return &transport_porosity;
        case Variable::vapour_pressure:
            return &vapour_pressure;
        case Variable::volumetric_mechanical_strain:
            return &volumetric_mechanical_strain;
        case Variable::volumetric_strain:
            return &volumetric_strain;
        default:
            OGS_FATAL(
                "No conversion to VariableType is provided for variable "
                "{:d}",
                static_cast<int>(v));
    };
}

static VariableArray::VariablePointer dropConst(
    VariableArray::VariablePointerConst const const_pointer)
{
    return std::visit(
        []<typename T>(T const* ptr) -> VariableArray::VariablePointer
        { return const_cast<T*>(ptr); },
        const_pointer);
}

VariableArray::VariablePointer VariableArray::address_of(Variable const v)
{
    return dropConst(const_cast<const VariableArray&>(*this).address_of(v));
}

bool VariableArray::is2D() const
{
    return maybeHasSize<5>(deformation_gradient) &&  //
           maybeHasSize<4>(mechanical_strain) &&     //
           maybeHasSize<4>(stress) &&                //
           maybeHasSize<4>(total_strain) &&          //
           maybeHasSize<4>(total_stress);
}

bool VariableArray::is3D() const
{
    return maybeHasSize<9>(deformation_gradient) &&  //
           maybeHasSize<6>(mechanical_strain) &&     //
           maybeHasSize<6>(stress) &&                //
           maybeHasSize<6>(total_strain) &&          //
           maybeHasSize<6>(total_stress);
}

}  // namespace MaterialPropertyLib
