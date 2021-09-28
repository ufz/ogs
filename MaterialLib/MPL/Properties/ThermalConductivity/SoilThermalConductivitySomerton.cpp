/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on February 17, 2021, 3:27 PM
 */

#include "SoilThermalConductivitySomerton.h"

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
//
// For 1D problems:
//
template <>
SoilThermalConductivitySomerton<1>::SoilThermalConductivitySomerton(
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

    double const lambda_try = dry_thermal_conductivity_(t, pos)[0];
    double const lambda_wet = wet_thermal_conductivity_(t, pos)[0];

    if (lambda_try > lambda_wet)
    {
        OGS_FATAL(
            "In 'SoilThermalConductivitySomerton', "
            "dry_thermal_conductivity of '{:g}' is larger than "
            "wet_thermal_conductivity of '{:g}'.",
            lambda_try, lambda_wet);
    }
}
template <>
PropertyDataType SoilThermalConductivitySomerton<1>::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& pos, double const t,
    double const /*dt*/) const
{
    double const S_L = std::get<double>(
        variable_array[static_cast<int>(Variable::liquid_saturation)]);

    if (S_L <= 0.0)
    {
        return dry_thermal_conductivity_(t, pos)[0];
    }

    if (S_L > 1.0)
    {
        return wet_thermal_conductivity_(t, pos)[0];
    }

    double const lambda_dry = dry_thermal_conductivity_(t, pos)[0];

    return lambda_dry +
           std::sqrt(S_L) * (wet_thermal_conductivity_(t, pos)[0] - lambda_dry);
}

template <>
PropertyDataType SoilThermalConductivitySomerton<1>::dValue(
    VariableArray const& variable_array, Variable const variable,
    ParameterLib::SpatialPosition const& pos, double const t,
    double const /*dt*/) const
{
    if (variable != Variable::liquid_saturation)
    {
        OGS_FATAL(
            "SoilThermalConductivitySomerton::dValue is implemented for "
            "derivatives with respect to liquid saturation only.");
    }

    double const S_L = std::get<double>(
        variable_array[static_cast<int>(Variable::liquid_saturation)]);

    if (S_L <= 0.0 || S_L > 1.0)
    {
        return 0.0;
    }

    double const lambda_dry = dry_thermal_conductivity_(t, pos)[0];
    double const lambda_wet = wet_thermal_conductivity_(t, pos)[0];

    return 0.5 * (lambda_wet - lambda_dry) / std::sqrt(S_L);
}

//
// For 2D and 3D problems:
//
template <int GlobalDimension>
SoilThermalConductivitySomerton<GlobalDimension>::
    SoilThermalConductivitySomerton(
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

    auto const lambda_try = dry_thermal_conductivity_(t, pos);
    auto const lambda_wet = wet_thermal_conductivity_(t, pos);

    if (lambda_try.size() != lambda_wet.size())
    {
        OGS_FATAL(
            "In 'SoilThermalConductivitySomerton' input data, the data size of "
            "dry_thermal_conductivity of '{:d}' is different from that of "
            "dry_thermal_conductivity of '{:d}'.",
            lambda_try.size(), lambda_wet.size());
    }

    for (std::size_t i = 0; i < lambda_try.size(); i++)
    {
        if (lambda_try[i] > lambda_wet[i])
        {
            OGS_FATAL(
                "In 'SoilThermalConductivitySomerton', "
                "dry_thermal_conductivity of '{:g}' is larger than "
                "wet_thermal_conductivity of '{:g}'.",
                lambda_try[i], lambda_wet[i]);
        }
    }
}

template <int GlobalDimension>
PropertyDataType SoilThermalConductivitySomerton<GlobalDimension>::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& pos, double const t,
    double const /*dt*/) const
{
    double const S_L = std::get<double>(
        variable_array[static_cast<int>(Variable::liquid_saturation)]);

    // (S_L <= 0.0)
    std::vector<double> lambda_data = dry_thermal_conductivity_(t, pos);

    if (S_L > 1.0)
    {
        lambda_data = wet_thermal_conductivity_(t, pos);
    }

    if (S_L > 0.0 && S_L <= 1.0)
    {
        auto const lambda_wet_data = wet_thermal_conductivity_(t, pos);

        for (std::size_t i = 0; i < lambda_data.size(); i++)
        {
            lambda_data[i] +=
                std::sqrt(S_L) * (lambda_wet_data[i] - lambda_data[i]);
        }
    }

    // Local coordinate transformation is only applied for the case that the
    // initial intrinsic permeability is given with orthotropic assumption.
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

template <int GlobalDimension>
PropertyDataType SoilThermalConductivitySomerton<GlobalDimension>::dValue(
    VariableArray const& variable_array, Variable const variable,
    ParameterLib::SpatialPosition const& pos, double const t,
    double const /*dt*/) const
{
    if (variable != Variable::liquid_saturation)
    {
        OGS_FATAL(
            "SoilThermalConductivitySomerton::dValue is implemented for "
            "derivatives with respect to liquid saturation only.");
    }

    double const S_L = std::get<double>(
        variable_array[static_cast<int>(Variable::liquid_saturation)]);

    if (S_L <= 0.0 || S_L > 1.0)
    {
        Eigen::Matrix<double, GlobalDimension, GlobalDimension> zero =
            Eigen::Matrix<double, GlobalDimension, GlobalDimension>::Zero();
        return zero;
    }

    auto const lambda_dry_data = dry_thermal_conductivity_(t, pos);
    auto const lambda_wet_data = wet_thermal_conductivity_(t, pos);

    std::vector<double> derivative_data;
    derivative_data.reserve(lambda_dry_data.size());
    for (std::size_t i = 0; i < lambda_dry_data.size(); i++)
    {
        derivative_data.emplace_back(
            0.5 * (lambda_wet_data[i] - lambda_dry_data[i]) / std::sqrt(S_L));
    }

    // Local coordinate transformation is only applied for the case that the
    // initial intrinsic permeability is given with orthotropic assumption.
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

template class SoilThermalConductivitySomerton<2>;
template class SoilThermalConductivitySomerton<3>;

}  // namespace MaterialPropertyLib
// namespace MaterialPropertyLib
