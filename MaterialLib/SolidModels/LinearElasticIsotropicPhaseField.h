/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

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
template <int DisplacementDim>
std::tuple<MathLib::KelvinVector::KelvinVectorType<
               DisplacementDim> /* sigma_real */,
           MathLib::KelvinVector::KelvinVectorType<
               DisplacementDim> /* sigma_tensile */,
           MathLib::KelvinVector::KelvinMatrixType<
               DisplacementDim> /* C_tensile */,
           MathLib::KelvinVector::KelvinMatrixType<
               DisplacementDim> /* C_compressive */,
           double /* strain_energy_tensile */,
           double /* elastic_energy */
           >
calculateDegradedStress(
    double const degradation,
    double const bulk_modulus,
    double const mu,
    MathLib::KelvinVector::KelvinVectorType<DisplacementDim> const& eps)
{
    static int const KelvinVectorSize =
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value;
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

    if (eps_curr_trace >= 0)
    {
        double const strain_energy_tensile =
            bulk_modulus / 2 * eps_curr_trace * eps_curr_trace +
            mu * epsd_curr.transpose() * epsd_curr;
        KelvinVector const sigma_tensile =
            bulk_modulus * eps_curr_trace * Invariants::identity2 +
            2 * mu * epsd_curr;
        C_tensile.template topLeftCorner<3, 3>().setConstant(bulk_modulus);
        C_tensile.noalias() += 2 * mu * P_dev * KelvinMatrix::Identity();
        double const elastic_energy = degradation * strain_energy_tensile;
        KelvinVector const sigma_real =
            degradation * sigma_tensile;  // + sigma_compressive, which is zero;
        return std::make_tuple(sigma_real, sigma_tensile, C_tensile,
                               C_compressive, strain_energy_tensile,
                               elastic_energy);
    }
    double const strain_energy_tensile = mu * epsd_curr.transpose() * epsd_curr;
    KelvinVector const sigma_tensile = 2 * mu * epsd_curr;
    KelvinVector const sigma_compressive =
        bulk_modulus * eps_curr_trace * Invariants::identity2;
    C_tensile.noalias() = 2 * mu * P_dev * KelvinMatrix::Identity();
    C_compressive.template topLeftCorner<3, 3>().setConstant(bulk_modulus);
    double const elastic_energy =
        bulk_modulus / 2 * eps_curr_trace * eps_curr_trace +
        mu * epsd_curr.transpose() * epsd_curr;
    KelvinVector const sigma_real =
        degradation * sigma_tensile + sigma_compressive;
    return std::make_tuple(sigma_real, sigma_tensile, C_tensile, C_compressive,
                           strain_energy_tensile, elastic_energy);
}
}  // namespace Phasefield
}  // namespace Solids
}  // namespace MaterialLib
