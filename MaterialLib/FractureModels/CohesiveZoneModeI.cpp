/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "CohesiveZoneModeI.h"
#include "LogPenalty.h"

#include "MathLib/MathTools.h"

namespace MaterialLib
{
namespace Fracture
{
namespace CohesiveZoneModeI
{
namespace
{
double computeDamage(double const damage_prev,
                     double const w_n,
                     double const w_np,
                     double const w_nf)
{
    return std::min(
        1.0,
        std::max(damage_prev, std::max(0.0, (w_n - w_np)) / (w_nf - w_np)));
}

}  // namespace

template <int DisplacementDim>
void CohesiveZoneModeI<DisplacementDim>::computeConstitutiveRelation(
    double const t,
    ParameterLib::SpatialPosition const& x,
    double const aperture0,
    Eigen::Ref<Eigen::VectorXd const>
        sigma0,
    Eigen::Ref<Eigen::VectorXd const>
    /*w_prev*/,
    Eigen::Ref<Eigen::VectorXd const>
        w,
    Eigen::Ref<Eigen::VectorXd const>
    /*sigma_prev*/,
    Eigen::Ref<Eigen::VectorXd>
        sigma,
    Eigen::Ref<Eigen::MatrixXd>
        C,
    typename FractureModelBase<DisplacementDim>::MaterialStateVariables&
        material_state_variables)
{
    assert(dynamic_cast<StateVariables<DisplacementDim> const*>(
               &material_state_variables) != nullptr);

    auto& state =
        static_cast<StateVariables<DisplacementDim>&>(material_state_variables);
    //reset damage in each iteration
    state.setInitialConditions();

    auto const mp = evaluatedMaterialProperties(t, x);

    C.setZero();

    // Separately compute shear and normal stresses because of the penalty for
    // the normal component.
    const int index_ns = DisplacementDim - 1;
    double const w_n = w[index_ns];
    for (int i = 0; i < index_ns; i++)
    {
        C(i, i) = mp.Ks;
    }

    sigma.noalias() = C * w;

    double const aperture = w_n + aperture0;

    sigma.coeffRef(index_ns) =
        mp.Kn * w_n * logPenalty(aperture0, aperture, _penalty_aperture_cutoff);

    C(index_ns, index_ns) =
        mp.Kn *
        logPenaltyDerivative(aperture0, aperture, _penalty_aperture_cutoff);

    sigma.noalias() += sigma0;

    // Exit if fracture is closing
    if (sigma[index_ns] < 0)
    {
        return;  /// Undamaged stiffness used in compression.
    }

    //
    // Continue with fracture opening.
    //
    // effective opening calculated by shift to state corresponding to sigma0=0
    const double w_n_effective = w_n + sigma0[index_ns] / mp.Kn;
    double degradation = 0.0;
    if (mp.w_nf == 0.0)  // can be used to model tension-cutoff
    {
        state.damage_prev = 1.0;
        state.damage = 1.0;
    }
    else
    {
        state.damage =
            computeDamage(state.damage_prev, w_n_effective, mp.w_np, mp.w_nf);
        degradation = ((1 - state.damage) * mp.w_np) /
                                   (mp.w_np + state.damage * (mp.w_nf - mp.w_np));
    }

    // Degrade stiffness tensor in tension.
    C = C * degradation;
    // Create effective relative displacement vector
    Eigen::VectorXd w_eff = w;
    // Only origin of normal entry is shifted, shear unaffected
    w_eff[index_ns] = w_n_effective;
    sigma = C * w_eff;

    if (state.damage > state.damage_prev)
    {
        // If damage is increasing, provide extension to consistent tangent.

        Eigen::Matrix<double, DisplacementDim, 1> dd_dw =
            Eigen::Matrix<double, DisplacementDim, 1>::Zero();

        double const tmp = mp.w_np + state.damage * (mp.w_nf - mp.w_np);
        dd_dw[index_ns] =
            (mp.w_np * mp.w_nf) / ((mp.w_nf - mp.w_np) * (tmp * tmp));

        C -= C * w_eff * (dd_dw).transpose();
    }

    // TODO (nagel) Initial stress not considered, yet.
    // sigma.noalias() += sigma0;
}

template class CohesiveZoneModeI<2>;
template class CohesiveZoneModeI<3>;

}  // namespace CohesiveZoneModeI
}  // namespace Fracture
}  // namespace MaterialLib
