/**
 * \file
 * \author Norbert Grunwald
 * \date   Jan, 2021
 * \brief
 *
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MaterialLib/MPL/Properties/IdealGasLawBinaryMixture.h"

#include "MaterialLib/PhysicalConstant.h"

namespace MaterialPropertyLib
{
IdealGasLawBinaryMixture::IdealGasLawBinaryMixture(std::string name)
{
    name_ = std::move(name);
}

PropertyDataType IdealGasLawBinaryMixture::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    const double R = MaterialLib::PhysicalConstant::IdealGasConstant;
    const double pGR = std::get<double>(
        variable_array[static_cast<int>(Variable::phase_pressure)]);
    const double T = std::get<double>(
        variable_array[static_cast<int>(Variable::temperature)]);
    const double MG = std::get<double>(
        variable_array[static_cast<int>(Variable::molar_mass)]);

    return pGR * MG / R / T;
}

PropertyDataType IdealGasLawBinaryMixture::dValue(
    VariableArray const& variable_array, Variable const variable,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    const double R = MaterialLib::PhysicalConstant::IdealGasConstant;
    const double pGR = std::get<double>(
        variable_array[static_cast<int>(Variable::phase_pressure)]);
    const double T = std::get<double>(
        variable_array[static_cast<int>(Variable::temperature)]);
    const double MG = std::get<double>(
        variable_array[static_cast<int>(Variable::molar_mass)]);
    const double dMG_dx = std::get<double>(
        variable_array[static_cast<int>(Variable::molar_mass_derivative)]);

    switch (variable)
    {
        case Variable::phase_pressure:
            return MG / R / T + pGR / R / T * dMG_dx;
        case Variable::capillary_pressure:
            return pGR / R / T * dMG_dx;
        case Variable::temperature:
            return pGR / R / T / T * (T * dMG_dx - MG);
        default:
            OGS_FATAL(
                "IdealGasLawBinaryMixture::dValue is implemented only for "
                "derivatives "
                "w.r.t. gas_pressure, capillary_pressure, and temperature.");
            break;
    }
    return 0.;
}

PropertyDataType IdealGasLawBinaryMixture::d2Value(
    VariableArray const& /*variable_array*/,
    Variable const /*primary_variable1*/, Variable const /*primary_variable2*/,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    OGS_FATAL("IdealGasLawBinaryMixture::d2Value is not implemented.");

    return 0.;
}

}  // namespace MaterialPropertyLib
