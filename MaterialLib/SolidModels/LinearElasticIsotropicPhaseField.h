/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "PhaseFieldExtension.h"
#include "LinearElasticIsotropic.h"
#include "KelvinVector.h"

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
        ProcessLib::KelvinVectorDimensions<DisplacementDim>::value;
    using KelvinVector = ProcessLib::KelvinVectorType<DisplacementDim>;
    using KelvinMatrix = ProcessLib::KelvinMatrixType<DisplacementDim>;

    explicit LinearElasticIsotropicPhaseField(
        typename LinearElasticIsotropic<
            DisplacementDim>::MaterialProperties&& material_properties)
        : LinearElasticIsotropic<DisplacementDim>(std::move(material_properties))
    {
    }

    std::unique_ptr<
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables>
    createMaterialStateVariables() override
    {
        return LinearElasticIsotropic<
            DisplacementDim>::createMaterialStateVariables();
    }

    boost::optional<std::tuple<KelvinVector,
                               std::unique_ptr<typename MechanicsBase<
                                   DisplacementDim>::MaterialStateVariables>,
                               KelvinMatrix>>
    integrateStress(
        double const t,
        ProcessLib::SpatialPosition const& x,
        double const dt,
        KelvinVector const& eps_prev,
        KelvinVector const& eps,
        KelvinVector const& sigma_prev,
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables const&
            material_state_variables) override
    {
        return LinearElasticIsotropic<DisplacementDim>::
                        integrateStress(t,
                                        x,
                                        dt,
                                        eps_prev,
                                        eps,
                                        sigma_prev,
                                        material_state_variables);
    }

    bool calculateDegradedStress(double const t,
                         ProcessLib::SpatialPosition const& x,
                         KelvinVector const& eps,
                         double& strain_energy_tensile,
                         KelvinVector& sigma_tensile,
                         KelvinVector& sigma_compressive,
                         KelvinMatrix& C_tensile,
                         KelvinMatrix& C_compressive,
                         KelvinVector& sigma_real,
                         double const degradation) const override
    {
        using Invariants =
            MaterialLib::SolidModels::Invariants<KelvinVectorSize>;
        // calculation of deviatoric parts
        auto const& P_dev = Invariants::deviatoric_projection;
        KelvinVector const epsd_curr = P_dev * eps;

        // Hydrostatic part for the stress and the tangent.
        double const eps_curr_trace = Invariants::trace(eps);

        auto const& K =
            LinearElasticIsotropic<DisplacementDim>::_mp.bulk_modulus(t, x);
        auto const& mu =
            LinearElasticIsotropic<DisplacementDim>::_mp.mu(t, x);

        C_tensile = KelvinMatrix::Zero();
        C_compressive = KelvinMatrix::Zero();

        if (eps_curr_trace >= 0)
        {
            strain_energy_tensile =
                 K / 2 * eps_curr_trace * eps_curr_trace +
                 mu * epsd_curr.transpose() * epsd_curr;
            sigma_tensile.noalias() =
                 K * eps_curr_trace * Invariants::identity2 +
                 2 * mu * epsd_curr;
            sigma_compressive.noalias() = KelvinVector::Zero();
            C_tensile.template topLeftCorner<3, 3>().setConstant(K);
            C_tensile.noalias() += 2 * mu * P_dev * KelvinMatrix::Identity();
        }
        else
        {
            strain_energy_tensile = mu * epsd_curr.transpose() * epsd_curr;
            sigma_tensile.noalias() = 2 * mu * epsd_curr;
            sigma_compressive.noalias() =
                 K * eps_curr_trace * Invariants::identity2;
            C_tensile.noalias() = 2 * mu * P_dev * KelvinMatrix::Identity();
            C_compressive.template topLeftCorner<3, 3>().setConstant(K);
        }

        sigma_real.noalias() = degradation * sigma_tensile + sigma_compressive;
        return true;
    }
};

extern template class LinearElasticIsotropicPhaseField<2>;
extern template class LinearElasticIsotropicPhaseField<3>;

}  // namespace Solids
}  // namespace MaterialLib
