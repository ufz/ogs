/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on November 10, 2020, 8:49 AM
 */

#include "StrainDependentPermeability.h"

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
StrainDependentPermeability<DisplacementDim>::StrainDependentPermeability(
    std::string name, ParameterLib::Parameter<double> const& k0,
    double const b1, double const b2, double const b3,
    double const minimum_permeability, double const maximum_permeability,
    ParameterLib::CoordinateSystem const* const local_coordinate_system)
    : k0_(k0),
      b1_(b1),
      b2_(b2),
      b3_(b3),
      minimum_permeability_(minimum_permeability),
      maximum_permeability_(maximum_permeability),
      local_coordinate_system_(local_coordinate_system)
{
    name_ = std::move(name);
}

template <int DisplacementDim>
void StrainDependentPermeability<DisplacementDim>::checkScale() const
{
    if (!std::holds_alternative<Medium*>(scale_))
    {
        OGS_FATAL(
            "The property 'StrainDependentPermeability' is "
            "implemented on the 'medium' scale only.");
    }
}

template <int DisplacementDim>
PropertyDataType StrainDependentPermeability<DisplacementDim>::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& pos, double const t,
    double const /*dt*/) const
{
    double const e_vol = std::get<double>(
        variable_array[static_cast<int>(Variable::volumetric_strain)]);
    double const e_vol_pls = std::get<double>(
        variable_array[static_cast<int>(Variable::equivalent_plastic_strain)]);

    auto k_data = k0_(t, pos);

    double const ten_base_exponent = e_vol > 0.0 ? b3_ * e_vol : b2_ * e_vol;
    double const factor =
        std::pow(10.0, ten_base_exponent) * std::exp(b1_ * e_vol_pls);

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
PropertyDataType StrainDependentPermeability<DisplacementDim>::dValue(
    VariableArray const& /*variable_array*/, Variable const /*variable*/,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    OGS_FATAL(
        "The derivative of the intrinsic permeability of "
        "StrainDependentPermeability"
        "is not implemented.");
}
template class StrainDependentPermeability<2>;
template class StrainDependentPermeability<3>;
}  // namespace MaterialPropertyLib
