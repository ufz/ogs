/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on February 17, 2021, 3:27 PM
 */

#include "SaturationWeightedThermalConductivity.h"

#include <cmath>
#include <limits>
#include <vector>

#include "BaseLib/Error.h"
#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Utils/CheckVanGenuchtenExponentRange.h"
#include "MaterialLib/MPL/Utils/FormEigenTensor.h"
#include "MaterialLib/MPL/Utils/GetSymmetricTensor.h"
#include "MaterialLib/MPL/VariableType.h"
#include "MathLib/KelvinVector.h"
#include "MathLib/MathTools.h"
#include "ParameterLib/CoordinateSystem.h"
#include "ParameterLib/Parameter.h"
#include "ParameterLib/SpatialPosition.h"

namespace MaterialPropertyLib
{
template <MeanType MeanType, int GlobalDimension>
SaturationWeightedThermalConductivity<MeanType, GlobalDimension>::
    SaturationWeightedThermalConductivity(
        std::string name,
        ParameterLib::Parameter<double> const& dry_thermal_conductivity,
        ParameterLib::Parameter<double> const& wet_thermal_conductivity,
        ParameterLib::CoordinateSystem const* const local_coordinate_system)
    : dry_thermal_conductivity_(dry_thermal_conductivity),
      wet_thermal_conductivity_(wet_thermal_conductivity),
      local_coordinate_system_(local_coordinate_system)
{
    name_ = std::move(name);

    ParameterLib::SpatialPosition const pos;
    double const t = std::numeric_limits<double>::quiet_NaN();

    auto const lambda_dry = dry_thermal_conductivity_(t, pos);
    auto const lambda_wet = wet_thermal_conductivity_(t, pos);

    if (lambda_dry.size() != lambda_wet.size())
    {
        OGS_FATAL(
            "In 'SaturationWeightedThermalConductivity' input data, the data "
            "size of "
            "dry_thermal_conductivity of '{:d}' is different from that of "
            "dry_thermal_conductivity of '{:d}'.",
            lambda_dry.size(), lambda_wet.size());
    }

    for (std::size_t i = 0; i < lambda_dry.size(); i++)
    {
        if (lambda_dry[i] > lambda_wet[i])
        {
            OGS_FATAL(
                "In 'SaturationWeightedThermalConductivity', "
                "dry_thermal_conductivity of '{:g}' is larger than "
                "wet_thermal_conductivity of '{:g}'.",
                lambda_dry[i], lambda_wet[i]);
        }
    }
    if constexpr (MeanType == MeanType::GEOMETRIC)
    {
        if (lambda_dry.size() != 1 &&
                lambda_dry.size() != GlobalDimension)
        {
            OGS_FATAL(
                    "The saturation weighted geometric mean"
                    "is not implemented for arbitrary anisotropic thermal "
                    "conductivities and requires to be in diagonal shape.");
        }
    }
}

template <MeanType MeanType, int GlobalDimension>
PropertyDataType SaturationWeightedThermalConductivity<MeanType, GlobalDimension>::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& pos, double const t,
    double const /*dt*/) const
{
    double const S_L = variable_array.liquid_saturation;

    // (S_L <= 0.0)
    std::vector<double> lambda_data = dry_thermal_conductivity_(t, pos);

    if (S_L >= 1.0)
    {
        lambda_data = wet_thermal_conductivity_(t, pos);
    }

    if (S_L > 0.0 && S_L <= 1.0)
    {

        for (std::size_t i = 0; i < lambda_data.size(); i++)
        {
                lambda_data[i] = getValue<MeanType>(S_L, lambda_data[i],
                                      wet_thermal_conductivity_(t, pos)[i]);
        }
    }

    if constexpr (GlobalDimension == 1)
    {
        return fromVector(lambda_data);
    }
    else
    {
        // Local coordinate transformation is only applied for the case that the
        // thermal conductivity is given with orthotropic assumption.
        if (local_coordinate_system_ && (lambda_data.size() == GlobalDimension))
        {
            Eigen::Matrix<double, GlobalDimension, GlobalDimension> const e =
                local_coordinate_system_->transformation<GlobalDimension>(pos);
            Eigen::Matrix<double, GlobalDimension, GlobalDimension> k =
            Eigen::Matrix<double, GlobalDimension, GlobalDimension>::Zero();

            for (int i = 0; i < GlobalDimension; ++i)
            {
                Eigen::Matrix<double, GlobalDimension, GlobalDimension> const
                    ei_otimes_ei = e.col(i) * e.col(i).transpose();

                k += lambda_data[i] * ei_otimes_ei;
            }
            return k;
        }

        return fromVector(lambda_data);
    }
}

template <MeanType MeanType, int GlobalDimension>
PropertyDataType SaturationWeightedThermalConductivity<MeanType, GlobalDimension>::dValue(
    VariableArray const& variable_array, Variable const variable,
    ParameterLib::SpatialPosition const& pos, double const t,
    double const /*dt*/) const
{
    if (variable != Variable::liquid_saturation)
    {
        OGS_FATAL(
            "SaturationWeightedThermalConductivity::dValue is implemented for "
            "derivatives with respect to liquid saturation only.");
    }

    double const S_L = variable_array.liquid_saturation;

    auto const lambda_dry_data = dry_thermal_conductivity_(t, pos);

    std::vector<double> derivative_data(lambda_dry_data.size(), 0.0);
    if (S_L <= 0.0 || S_L > 1.0)
    {
        return fromVector(derivative_data);
    }
    for (std::size_t i = 0; i < lambda_dry_data.size(); i++)
    {
            derivative_data[i] = getDValue<MeanType>(S_L, lambda_dry_data[i],
                                      wet_thermal_conductivity_(t, pos)[i]);
    }
    if constexpr (GlobalDimension == 1)
    {
        return fromVector(derivative_data);
    }
    else
    {
        // Local coordinate transformation is only applied for the case that the
        // thermal conductivity is given with orthotropic assumption.
        if (local_coordinate_system_ && (derivative_data.size() == GlobalDimension))
        {
            Eigen::Matrix<double, GlobalDimension, GlobalDimension> const e =
                local_coordinate_system_->transformation<GlobalDimension>(pos);
            Eigen::Matrix<double, GlobalDimension, GlobalDimension> k =
            Eigen::Matrix<double, GlobalDimension, GlobalDimension>::Zero();

            for (int i = 0; i < GlobalDimension; ++i)
            {
                Eigen::Matrix<double, GlobalDimension, GlobalDimension> const
                    ei_otimes_ei = e.col(i) * e.col(i).transpose();

                k += derivative_data[i] * ei_otimes_ei;
            }
            return k;
        }

        return fromVector(derivative_data);
    }
}

template class SaturationWeightedThermalConductivity<MeanType::ARITHMETIC_LINEAR, 1>;
template class SaturationWeightedThermalConductivity<MeanType::ARITHMETIC_SQUAREROOT, 1>;
template class SaturationWeightedThermalConductivity<MeanType::GEOMETRIC, 1>;
template class SaturationWeightedThermalConductivity<MeanType::ARITHMETIC_LINEAR, 2>;
template class SaturationWeightedThermalConductivity<MeanType::ARITHMETIC_SQUAREROOT, 2>;
template class SaturationWeightedThermalConductivity<MeanType::GEOMETRIC, 2>;
template class SaturationWeightedThermalConductivity<MeanType::ARITHMETIC_LINEAR, 3>;
template class SaturationWeightedThermalConductivity<MeanType::ARITHMETIC_SQUAREROOT, 3>;
template class SaturationWeightedThermalConductivity<MeanType::GEOMETRIC, 3>;

}  // namespace MaterialPropertyLib
// namespace MaterialPropertyLib
