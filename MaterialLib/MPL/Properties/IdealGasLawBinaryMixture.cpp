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

// Enums with indices for the sorting of the partial density and its derivatives
// in the result vector or the result tensor of the functions value and dValue.
enum class VariableIndex : int
{
    gas_pressure,
    capillary_pressure,
    temperature
};
enum class DensityIndex : int
{
    gas_phase_density,
    vapour_density,
    air_density
};

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
    const double MW = std::get<double>(
        variable_array[static_cast<int>(Variable::molar_mass_vapour)]);
    const double MC = std::get<double>(
        variable_array[static_cast<int>(Variable::molar_mass)]);
    // mole fraction of the vapour component
    const double xnWG = std::get<double>(
        variable_array[static_cast<int>(Variable::molar_fraction)]);
    const double xnCG = 1. - xnWG;

    const double rhoWGR = xnWG * pGR * MW / R / T;
    const double rhoCGR = xnCG * pGR * MC / R / T;
    Eigen::Matrix<double, 3, 1> result;

    result[static_cast<int>(DensityIndex::vapour_density)] = rhoWGR;
    result[static_cast<int>(DensityIndex::air_density)] = rhoCGR;
    result[static_cast<int>(DensityIndex::gas_phase_density)] = rhoWGR + rhoCGR;

    return result;
}

PropertyDataType IdealGasLawBinaryMixture::dValue(
    VariableArray const& variable_array, Variable const /*primary_variable*/,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    const double R = MaterialLib::PhysicalConstant::IdealGasConstant;
    const double pGR = std::get<double>(
        variable_array[static_cast<int>(Variable::phase_pressure)]);
    const double T = std::get<double>(
        variable_array[static_cast<int>(Variable::temperature)]);
    const double MW = std::get<double>(
        variable_array[static_cast<int>(Variable::molar_mass_vapour)]);
    const double MC = std::get<double>(
        variable_array[static_cast<int>(Variable::molar_mass)]);
    // mole fraction of the vapour component
    const double xnWG = std::get<double>(
        variable_array[static_cast<int>(Variable::molar_fraction)]);
    const double xnCG = 1. - xnWG;
    const double dxnWG_dpGR =
        std::get<double>(variable_array[static_cast<int>(Variable::dxn_dpGR)]);
    const double dxnWG_dpCap =
        std::get<double>(variable_array[static_cast<int>(Variable::dxn_dpCap)]);
    const double dxnWG_dT =
        std::get<double>(variable_array[static_cast<int>(Variable::dxn_dT)]);

    Eigen::Matrix<double, 3, 3> result;

    const double drhoCGR_dpGR = -MC / R / T * (pGR * dxnWG_dpGR - xnCG);
    const double drhoCGR_dpCap = -pGR * MC / R / T * dxnWG_dpCap;
    const double drhoCGR_dT = -pGR * MC / R / T * (dxnWG_dT + xnCG / T);

    const double drhoWGR_dpGR = MW / R / T * (pGR * dxnWG_dpGR + xnWG);
    const double drhoWGR_dpCap = pGR * MW / R / T * dxnWG_dpCap;
    const double drhoWGR_dT = pGR * MW / R / T / T * (T * dxnWG_dT - xnWG);

    const double drhoGR_dpGR = drhoCGR_dpGR + drhoWGR_dpGR;
    const double drhoGR_dpCap = drhoCGR_dpCap + drhoWGR_dpCap;
    const double drhoGR_dT = drhoCGR_dT + drhoWGR_dT;

    result(static_cast<int>(DensityIndex::air_density),
           static_cast<int>(VariableIndex::gas_pressure)) = drhoCGR_dpGR;
    result(static_cast<int>(DensityIndex::air_density),
           static_cast<int>(VariableIndex::capillary_pressure)) = drhoCGR_dpCap;
    result(static_cast<int>(DensityIndex::air_density),
           static_cast<int>(VariableIndex::temperature)) = drhoCGR_dT;

    result(static_cast<int>(DensityIndex::vapour_density),
           static_cast<int>(VariableIndex::gas_pressure)) = drhoWGR_dpGR;
    result(static_cast<int>(DensityIndex::vapour_density),
           static_cast<int>(VariableIndex::capillary_pressure)) = drhoWGR_dpCap;
    result(static_cast<int>(DensityIndex::vapour_density),
           static_cast<int>(VariableIndex::temperature)) = drhoWGR_dT;

    result(static_cast<int>(DensityIndex::gas_phase_density),
           static_cast<int>(VariableIndex::gas_pressure)) = drhoGR_dpGR;
    result(static_cast<int>(DensityIndex::gas_phase_density),
           static_cast<int>(VariableIndex::capillary_pressure)) = drhoGR_dpCap;
    result(static_cast<int>(DensityIndex::gas_phase_density),
           static_cast<int>(VariableIndex::temperature)) = drhoGR_dT;

    return result;
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
