/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <spdlog/spdlog.h>

#include <Eigen/Eigenvalues>

#include "MaterialLib/SolidModels/Ehlers.h"

namespace ProcessLib
{
namespace SmallDeformationNonlocal
{
/// Computes the damage internal material variable explicitly based on the
/// results obtained from the local stress return algorithm.
inline double calculateDamage(double const kappa_d, double const alpha_d,
                              double const beta_d)
{
    double const damage = (1 - beta_d) * (1 - std::exp(-kappa_d / alpha_d));

    if (damage < 0. || damage >= 1.)
    {
        OGS_FATAL("Damage value {:g} outside of [0,1) interval.", damage);
    }

    return damage;
}

template <int DisplacementDim, typename KelvinVectorType>
double calculateDamageKappaD(
    double const eps_p_eff_diff,
    KelvinVectorType const& sigma,
    double const kappa_d_prev,
    double const h_d,
    MaterialLib::Solids::Ehlers::MaterialProperties const& mp)
{
    // Default case of the rate problem. Updated below if volumetric plastic
    // strain rate is positive (dilatancy).

    // non-const for Eigen solver.
    auto sigma_tensor = MathLib::KelvinVector::kelvinVectorToTensor(sigma);

    Eigen::EigenSolver<decltype(sigma_tensor)> eigen_solver(sigma_tensor);
    auto const principal_stress = real(eigen_solver.eigenvalues().array());
    double const prod_stress = std::sqrt(principal_stress.square().sum());

    // Brittleness decrease with confinement for the nonlinear flow rule.
    // ATTENTION: For linear flow rule -> constant brittleness.
    double const tensile_strength =
        std::sqrt(3.0) * mp.kappa / (1 + std::sqrt(3.0) * mp.beta);
    double const r_s = prod_stress / tensile_strength - 1.;

    // Compute normalizing strain.
    double const x_s = [](double const h_d, double const r_s) {
        if (r_s < 0)
        {
            return 1.;
        }
        if (r_s <= 1)
        {
            return 1. + h_d * r_s * r_s;
        }
        return 1. - 3 * h_d + 4 * h_d * std::sqrt(r_s);
    }(h_d, r_s);

    return kappa_d_prev + eps_p_eff_diff / x_s;
}
}  // namespace SmallDeformationNonlocal
}  // namespace ProcessLib
