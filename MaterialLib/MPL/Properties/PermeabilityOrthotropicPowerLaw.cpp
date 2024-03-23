/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "MaterialLib/MPL/Properties/PermeabilityOrthotropicPowerLaw.h"

#include <algorithm>
#include <cmath>

#include "MaterialLib/MPL/Medium.h"
#include "MathLib/KelvinVector.h"
#include "ParameterLib/CoordinateSystem.h"

namespace MaterialPropertyLib
{
template <int DisplacementDim>
PermeabilityOrthotropicPowerLaw<DisplacementDim>::
    PermeabilityOrthotropicPowerLaw(
        std::string name,
        std::array<double, DisplacementDim>
            intrinsic_permeabilities,
        std::array<double, DisplacementDim>
            exponents,
        ParameterLib::CoordinateSystem const* const local_coordinate_system)
    : k_(std::move(intrinsic_permeabilities)),
      lambda_(std::move(exponents)),
      local_coordinate_system_(local_coordinate_system)
{
    name_ = std::move(name);
}

template <int DisplacementDim>
void PermeabilityOrthotropicPowerLaw<DisplacementDim>::checkScale() const
{
    if (!std::holds_alternative<Medium*>(scale_))
    {
        OGS_FATAL(
            "The property 'PermeabilityOrthotropicPowerLaw' is implemented on "
            "the 'medium' scales only.");
    }
}
template <int DisplacementDim>
PropertyDataType PermeabilityOrthotropicPowerLaw<DisplacementDim>::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& pos, double const /*t*/,
    double const /*dt*/) const
{
    auto const phi = variable_array.transport_porosity;
    // TODO (naumov) The phi0 must be evaluated once upon
    // creation/initialization and be stored in a local state.
    // For now assume porosity's initial value does not change with time.
    auto const medium = std::get<Medium*>(scale_);
    auto const phi_0 =
        medium->hasProperty(PropertyType::transport_porosity)
            ? medium->property(PropertyType::transport_porosity)
                  .template initialValue<double>(
                      pos, std::numeric_limits<double>::quiet_NaN())
            : medium->property(PropertyType::porosity)
                  .template initialValue<double>(
                      pos, std::numeric_limits<double>::quiet_NaN());

    Eigen::Matrix<double, DisplacementDim, DisplacementDim> k =
        Eigen::Matrix<double, DisplacementDim, DisplacementDim>::Zero();

    Eigen::Matrix<double, DisplacementDim, DisplacementDim> const e =
        local_coordinate_system_ == nullptr
            ? Eigen::Matrix<double, DisplacementDim,
                            DisplacementDim>::Identity()
            : local_coordinate_system_->transformation<DisplacementDim>(pos);

    // k = \sum_i k_i (\phi / \phi_0)^{\lambda_i} e_i \otimes e_i
    // e_i \otimes e_i = square matrix e_i,0^2 e_i,0*e_i,1 etc.
    for (int i = 0; i < DisplacementDim; ++i)
    {
        Eigen::Matrix<double, DisplacementDim, DisplacementDim> const
            ei_otimes_ei = e.col(i) * e.col(i).transpose();

        k += k_[i] * std::pow(phi / phi_0, lambda_[i]) * ei_otimes_ei;
    }
    return k;
}
template <int DisplacementDim>
PropertyDataType PermeabilityOrthotropicPowerLaw<DisplacementDim>::dValue(
    VariableArray const& /*variable_array*/, Variable const /*variable*/,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    OGS_FATAL("PermeabilityOrthotropicPowerLaw::dValue is not implemented.");
}

template class PermeabilityOrthotropicPowerLaw<2>;
template class PermeabilityOrthotropicPowerLaw<3>;
}  // namespace MaterialPropertyLib
