/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on June 4, 2020, 10:13 AM
 */

#include "PermeabilityMohrCoulombFailureIndexModel.h"

#include <boost/math/constants/constants.hpp>
#include <cmath>
#include <limits>

#include "BaseLib/Error.h"
#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Utils/FormEigenTensor.h"
#include "MaterialLib/MPL/Utils/GetSymmetricTensor.h"
#include "MathLib/KelvinVector.h"
#include "MathLib/MathTools.h"
#include "ParameterLib/CoordinateSystem.h"
#include "ParameterLib/Parameter.h"

namespace MaterialPropertyLib
{
template <int DisplacementDim>
PermeabilityMohrCoulombFailureIndexModel<DisplacementDim>::
    PermeabilityMohrCoulombFailureIndexModel(
        std::string name, ParameterLib::Parameter<double> const& k0,
        double const kr, double const b, double const c, double const phi,
        double const k_max, double const t_sigma_max,
        ParameterLib::CoordinateSystem const* const local_coordinate_system)
    : k0_(k0),
      kr_(kr),
      b_(b),
      c_(c),
      phi_(boost::math::constants::degree<double>() * phi),
      k_max_(k_max),
      t_sigma_max_(t_sigma_max),
      local_coordinate_system_(local_coordinate_system)
{
    const double t_sigma_upper = c_ / std::tan(phi_);
    if (t_sigma_max_ <= 0.0 || t_sigma_max_ > t_sigma_upper ||
        std::fabs(t_sigma_max_ - t_sigma_upper) <
            std::numeric_limits<double>::epsilon())
    {
        OGS_FATAL(
            "Tensile strength parameter of {:e} is out of the range (0, "
            "c/tan(phi)) = (0, {:e})",
            t_sigma_max_, t_sigma_upper);
    }

    name_ = std::move(name);
}

template <int DisplacementDim>
void PermeabilityMohrCoulombFailureIndexModel<DisplacementDim>::checkScale()
    const
{
    if (!std::holds_alternative<Medium*>(scale_))
    {
        OGS_FATAL(
            "The property 'PermeabilityMohrCoulombFailureIndexModel' is "
            "implemented on the 'medium' scale only.");
    }
}

template <int DisplacementDim>
PropertyDataType
PermeabilityMohrCoulombFailureIndexModel<DisplacementDim>::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& pos, double const t,
    double const /*dt*/) const
{
    auto const& stress_vector = std::get<SymmetricTensor<DisplacementDim>>(
        variable_array[static_cast<int>(Variable::total_stress)]);

    auto const& stress_tensor =
        formEigenTensor<3>(static_cast<PropertyDataType>(stress_vector));

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 3, 3>>
        eigenvalue_solver(stress_tensor);

    // Principle stress
    auto const sigma = eigenvalue_solver.eigenvalues();

    auto k_data = k0_(t, pos);

    double const max_sigma = std::max(std::fabs(sigma[0]), std::fabs(sigma[2]));

    if (max_sigma < std::numeric_limits<double>::epsilon())
    {
        return fromVector(k_data);
    }

    double const sigma_m = 0.5 * (sigma[2] + sigma[0]);

    double const tau_m = 0.5 * std::fabs(sigma[2] - sigma[0]);
    double f = 0.0;
    if (sigma_m > t_sigma_max_)
    {
        // tensile failure criterion
        f = sigma_m / t_sigma_max_;

        double const tau_tt =
            c_ * std::cos(phi_) - t_sigma_max_ * std::sin(phi_);

        f = std::max(f, tau_m / tau_tt);
    }
    else
    {
        // Mohr Coulomb failure criterion
        f = tau_m / (c_ * std::cos(phi_) - sigma_m * std::sin(phi_));
    }

    if (f >= 1.0)
    {
        const double exp_value = std::exp(b_ * f);
        for (auto& k_i : k_data)
        {
            k_i = std::min(k_i + kr_ * exp_value, k_max_);
        }
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
PropertyDataType
PermeabilityMohrCoulombFailureIndexModel<DisplacementDim>::dValue(
    VariableArray const& /*variable_array*/, Variable const /*variable*/,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    OGS_FATAL(
        "The derivative of the intrinsic permeability k(sigma, ...) with "
        "respect to stress tensor (sigma) is not implemented because that "
        "dk/du is normally omitted.");
}
template class PermeabilityMohrCoulombFailureIndexModel<2>;
template class PermeabilityMohrCoulombFailureIndexModel<3>;
}  // namespace MaterialPropertyLib
