/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on January 12, 2021, x:xx AM
 */

#include "GasPressureDependentPermeability.h"

#include <cmath>
#include <limits>

#include "BaseLib/Error.h"
#include "MaterialLib/MPL/Medium.h"
#include "MathLib/MathTools.h"
#include "ParameterLib/CoordinateSystem.h"
#include "ParameterLib/Parameter.h"

namespace MaterialPropertyLib
{
template <int DisplacementDim>
GasPressureDependentPermeability<DisplacementDim>::
    GasPressureDependentPermeability(
        std::string name, ParameterLib::Parameter<double> const& k0,
        double const a1, double const a2, double const pressure_threshold,
        double const minimum_permeability, double const maximum_permeability,
        ParameterLib::CoordinateSystem const* const local_coordinate_system)
    : k0_(k0),
      a1_(a1),
      a2_(a2),
      pressure_threshold_(pressure_threshold),
      minimum_permeability_(minimum_permeability),
      maximum_permeability_(maximum_permeability),
      local_coordinate_system_(local_coordinate_system)
{
    name_ = std::move(name);
}

template <int DisplacementDim>
void GasPressureDependentPermeability<DisplacementDim>::checkScale() const
{
    if (!std::holds_alternative<Medium*>(scale_))
    {
        OGS_FATAL(
            "The property 'GasPressureDependentPermeability' is implemented on "
            "the 'medium' scale only.");
    }
}

template <int DisplacementDim>
PropertyDataType GasPressureDependentPermeability<DisplacementDim>::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& pos, double const t,
    double const /*dt*/) const
{
    double const gas_pressure = std::get<double>(
        variable_array[static_cast<int>(Variable::phase_pressure)]);

    auto k_data = k0_(t, pos);

    double const factor =
        (gas_pressure <= pressure_threshold_)
            ? (1.0 + a1_ * gas_pressure / 1.0e6)
            : (a2_ * (gas_pressure - pressure_threshold_) / 1.0e6 + 1.0 +
               a1_ * pressure_threshold_ / 1.0e6);

    for (auto& k_i : k_data)
    {
        k_i = std::clamp(k_i * factor, minimum_permeability_,
                         maximum_permeability_);
    }

    // Local coordinate transformation is only applied for the case that the
    // initial intrinsic permeability is given with orthotropic assumption.
    if (local_coordinate_system_ && (k_data.size() == DisplacementDim))
    {
        Eigen::Matrix<double, DisplacementDim, DisplacementDim> const e =
            local_coordinate_system_->transformation<DisplacementDim>(pos);
        Eigen::Matrix<double, DisplacementDim, DisplacementDim> k =
            Eigen::Matrix<double, DisplacementDim, DisplacementDim>::Zero();

        for (int i = 0; i < DisplacementDim; ++i)
        {
            Eigen::Matrix<double, DisplacementDim, DisplacementDim> const
                ei_otimes_ei = e.col(i) * e.col(i).transpose();

            k += k_data[i] * ei_otimes_ei;
        }
        return k;
    }

    return fromVector(k_data);
}

template <int DisplacementDim>
PropertyDataType GasPressureDependentPermeability<DisplacementDim>::dValue(
    VariableArray const& /*variable_array*/, Variable const /*variable*/,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    OGS_FATAL(
        "The derivative of the intrinsic permeability of "
        "GasPressureDependentPermeability"
        "is not implemented.");
}

template class GasPressureDependentPermeability<2>;
template class GasPressureDependentPermeability<3>;
}  // namespace MaterialPropertyLib
