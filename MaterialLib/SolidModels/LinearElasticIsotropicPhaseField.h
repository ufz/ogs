/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MathLib/KelvinVector.h"

#include "LinearElasticIsotropic.h"
#include "PhaseFieldExtension.h"

namespace MaterialLib
{
namespace Solids
{
template <int DisplacementDim>
class LinearElasticIsotropicPhaseField final
    : public LinearElasticIsotropic<DisplacementDim>,
      public PhaseFieldExtension<DisplacementDim>
{
public:
    static int const KelvinVectorSize =
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value;
    using KelvinVector =
        MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;
    using KelvinMatrix =
        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>;

    explicit LinearElasticIsotropicPhaseField(
        typename LinearElasticIsotropic<DisplacementDim>::MaterialProperties&&
            material_properties)
        : LinearElasticIsotropic<DisplacementDim>(
              std::move(material_properties))
    {
    }

    std::unique_ptr<
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables>
    createMaterialStateVariables() const override
    {
        return LinearElasticIsotropic<
            DisplacementDim>::createMaterialStateVariables();
    }

    boost::optional<std::tuple<KelvinVector,
                               std::unique_ptr<typename MechanicsBase<
                                   DisplacementDim>::MaterialStateVariables>,
                               KelvinMatrix>>
    integrateStress(
        double const t, ProcessLib::SpatialPosition const& x, double const dt,
        KelvinVector const& eps_prev, KelvinVector const& eps,
        KelvinVector const& sigma_prev,
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables const&
            material_state_variables,
        double const T) const override
    {
        return LinearElasticIsotropic<DisplacementDim>::integrateStress(
            t, x, dt, eps_prev, eps, sigma_prev, material_state_variables, T);
    }

    double computeFreeEnergyDensity(
        double const t,
        ProcessLib::SpatialPosition const& x,
        double const dt,
        KelvinVector const& eps,
        KelvinVector const& sigma,
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables const&
            material_state_variables) const override
    {
        return LinearElasticIsotropic<DisplacementDim>::
            computeFreeEnergyDensity(
                t, x, dt, eps, sigma, material_state_variables);
    }

    /** Decompose the stiffness into tensile and compressive part.
     * Judging by the physical observations, compression perpendicular
     * to a crack does not cause crack propagation. Thus,
     * the phase-field parameter is only involved into the tensile part
     * to degrade the elastic strain energy.
     */
    bool calculateDegradedStress(double const t,
                                 ProcessLib::SpatialPosition const& x,
                                 KelvinVector const& eps,
                                 double& strain_energy_tensile,
                                 KelvinVector& sigma_tensile,
                                 KelvinVector& sigma_compressive,
                                 KelvinMatrix& C_tensile,
                                 KelvinMatrix& C_compressive,
                                 KelvinVector& sigma_real,
                                 double const degradation,
                                 double& elastic_energy) const override
    {
        using Invariants = MathLib::KelvinVector::Invariants<KelvinVectorSize>;
        // calculation of deviatoric parts
        auto const& P_dev = Invariants::deviatoric_projection;
        KelvinVector const epsd_curr = P_dev * eps;

        // Hydrostatic part for the stress and the tangent.
        double const eps_curr_trace = Invariants::trace(eps);

        auto const& K =
            LinearElasticIsotropic<DisplacementDim>::_mp.bulk_modulus(t, x);
        auto const& mu = LinearElasticIsotropic<DisplacementDim>::_mp.mu(t, x);

        C_tensile = KelvinMatrix::Zero();
        C_compressive = KelvinMatrix::Zero();

        if (eps_curr_trace >= 0)
        {
            strain_energy_tensile = K / 2 * eps_curr_trace * eps_curr_trace +
                                    mu * epsd_curr.transpose() * epsd_curr;
            sigma_tensile.noalias() =
                K * eps_curr_trace * Invariants::identity2 + 2 * mu * epsd_curr;
            sigma_compressive.noalias() = KelvinVector::Zero();
            C_tensile.template topLeftCorner<3, 3>().setConstant(K);
            C_tensile.noalias() += 2 * mu * P_dev * KelvinMatrix::Identity();
            elastic_energy = degradation * strain_energy_tensile;
        }
        else
        {
            strain_energy_tensile = mu * epsd_curr.transpose() * epsd_curr;
            sigma_tensile.noalias() = 2 * mu * epsd_curr;
            sigma_compressive.noalias() =
                K * eps_curr_trace * Invariants::identity2;
            C_tensile.noalias() = 2 * mu * P_dev * KelvinMatrix::Identity();
            C_compressive.template topLeftCorner<3, 3>().setConstant(K);
            elastic_energy = K / 2 * eps_curr_trace * eps_curr_trace +
                             mu * epsd_curr.transpose() * epsd_curr;
        }

        sigma_real.noalias() = degradation * sigma_tensile + sigma_compressive;
        return true;
    }
};

extern template class LinearElasticIsotropicPhaseField<2>;
extern template class LinearElasticIsotropicPhaseField<3>;

}  // namespace Solids
}  // namespace MaterialLib
