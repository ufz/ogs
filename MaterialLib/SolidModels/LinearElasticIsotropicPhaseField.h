/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <boost/math/special_functions/pow.hpp>

#include "MathLib/KelvinVector.h"

namespace MaterialLib
{
namespace Solids
{
namespace Phasefield
{
/** Decompose the stiffness into tensile and compressive part.
 * Judging by the physical observations, compression perpendicular
 * to a crack does not cause crack propagation. Thus,
 * the phase-field parameter is only involved into the tensile part
 * to degrade the elastic strain energy.
 */

/// heaviside function returns 1.0 if the argument is positive and 0.0 if
/// negative
inline double heaviside(double const v)
{
    return (v < 0) ? 0.0 : 1.0;
}

/// Macaulay brackets: positive strain is tensile and negative strain for
/// compressive
inline double macaulayTensile(double const v)
{
    return v * heaviside(v);
}
inline double macaulayCompressive(double v)
{
    return v * (1 - heaviside(v));
}

template <int DisplacementDim>
std::tuple<MathLib::KelvinVector::KelvinVectorType<
               DisplacementDim> /* sigma_real */,
           MathLib::KelvinVector::KelvinVectorType<
               DisplacementDim> /* sigma_tensile */,
           MathLib::KelvinVector::KelvinMatrixType<
               DisplacementDim> /* C_tensile */,
           MathLib::KelvinVector::KelvinMatrixType<
               DisplacementDim> /* C_compressive */,
           double /* strain_energy_tensile */, double /* elastic_energy */
           >
calculateVolDevDegradedStress(
    double const degradation, double const bulk_modulus, double const mu,
    MathLib::KelvinVector::KelvinVectorType<DisplacementDim> const& eps)
{
    static int const KelvinVectorSize =
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);
    using KelvinVector =
        MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;
    using KelvinMatrix =
        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>;
    using Invariants = MathLib::KelvinVector::Invariants<KelvinVectorSize>;
    // calculation of deviatoric parts
    auto const& P_dev = Invariants::deviatoric_projection;
    KelvinVector const epsd_curr = P_dev * eps;

    // Hydrostatic part for the stress and the tangent.
    double const eps_curr_trace = Invariants::trace(eps);

    KelvinMatrix C_tensile = KelvinMatrix::Zero();
    KelvinMatrix C_compressive = KelvinMatrix::Zero();

    auto strain_energy_computation_vol = [&](auto&& macaulay)
    {
        auto macaulay_squared = [&macaulay](double x)
        { return boost::math::pow<2>(macaulay(x)); };
        return bulk_modulus / 2 * macaulay_squared(eps_curr_trace);
    };

    auto stress_computation_vol = [&](auto&& macaulay)
    { return bulk_modulus * macaulay(eps_curr_trace) * Invariants::identity2; };

    auto hs = [&](double const v) { return heaviside(v); };

    auto mt = [&](double const v) { return macaulayTensile(v); };

    auto mc = [&](double const v) { return macaulayCompressive(v); };

    double const strain_energy_tensile = strain_energy_computation_vol(mt) +
                                         mu * epsd_curr.transpose() * epsd_curr;

    KelvinVector const sigma_tensile =
        stress_computation_vol(mt) + 2 * mu * epsd_curr;

    KelvinVector const sigma_compressive = stress_computation_vol(mc);

    C_tensile.template topLeftCorner<3, 3>().setConstant(bulk_modulus *
                                                         hs(eps_curr_trace));
    C_tensile.noalias() += 2 * mu * P_dev * KelvinMatrix::Identity();

    C_compressive.template topLeftCorner<3, 3>().setConstant(
        bulk_modulus * (1 - hs(eps_curr_trace)));

    double const elastic_energy =
        bulk_modulus / 2 * eps_curr_trace * eps_curr_trace +
        mu * epsd_curr.transpose() * epsd_curr;
    KelvinVector const sigma_real =
        degradation * sigma_tensile + sigma_compressive;
    return std::make_tuple(sigma_real, sigma_tensile, C_tensile, C_compressive,
                           strain_energy_tensile, elastic_energy);
}

template <int DisplacementDim>
std::tuple<MathLib::KelvinVector::KelvinVectorType<
               DisplacementDim> /* sigma_real */,
           MathLib::KelvinVector::KelvinVectorType<
               DisplacementDim> /* sigma_tensile */,
           MathLib::KelvinVector::KelvinMatrixType<
               DisplacementDim> /* C_tensile */,
           MathLib::KelvinVector::KelvinMatrixType<
               DisplacementDim> /* C_compressive */,
           double /* strain_energy_tensile */, double /* elastic_energy */
           >
calculateIsotropicDegradedStress(
    double const degradation,
    double const bulk_modulus,
    double const mu,
    MathLib::KelvinVector::KelvinVectorType<DisplacementDim> const& eps)
{
    static int const KelvinVectorSize =
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);
    using KelvinVector =
        MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;
    using KelvinMatrix =
        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>;
    using Invariants = MathLib::KelvinVector::Invariants<KelvinVectorSize>;
    // calculation of deviatoric parts
    auto const& P_dev = Invariants::deviatoric_projection;
    KelvinVector const epsd_curr = P_dev * eps;

    // Hydrostatic part for the stress and the tangent.
    double const eps_curr_trace = Invariants::trace(eps);

    KelvinMatrix C_tensile = KelvinMatrix::Zero();
    KelvinMatrix C_compressive = KelvinMatrix::Zero();

    double const strain_energy_tensile =
        bulk_modulus / 2 * eps_curr_trace * eps_curr_trace +
        mu * epsd_curr.transpose() * epsd_curr;
    double const elastic_energy = degradation * strain_energy_tensile;
    KelvinVector const sigma_tensile =
        bulk_modulus * eps_curr_trace * Invariants::identity2 +
        2 * mu * epsd_curr;
    C_tensile.template topLeftCorner<3, 3>().setConstant(bulk_modulus);
    C_tensile.noalias() += 2 * mu * P_dev * KelvinMatrix::Identity();
    KelvinVector const sigma_real = degradation * sigma_tensile;

    return std::make_tuple(sigma_real, sigma_tensile, C_tensile, C_compressive,
                           strain_energy_tensile, elastic_energy);
}

}  // namespace Phasefield
}  // namespace Solids
}  // namespace MaterialLib
