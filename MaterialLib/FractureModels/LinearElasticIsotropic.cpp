/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LinearElasticIsotropic.h"

#include "LogPenalty.h"

namespace MaterialLib
{
namespace Fracture
{
template <int DisplacementDim>
void LinearElasticIsotropic<DisplacementDim>::computeConstitutiveRelation(
    double const t,
    ProcessLib::SpatialPosition const& x,
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
    material_state_variables.reset();

    const int index_ns = DisplacementDim - 1;
    C.setZero();
    for (int i = 0; i < index_ns; i++)
        C(i, i) = _mp.shear_stiffness(t, x)[0];

    sigma.noalias() = C * w;

    double const aperture = w[index_ns] + aperture0;

    sigma.coeffRef(index_ns) =
        _mp.normal_stiffness(t, x)[0] * w[index_ns] *
        logPenalty(aperture0, aperture, _penalty_aperture_cutoff);

    C(index_ns, index_ns) =
        _mp.normal_stiffness(t, x)[0] *
        logPenaltyDerivative(aperture0, aperture, _penalty_aperture_cutoff);

    sigma.noalias() += sigma0;

    // correction for an opening fracture
    if (_tension_cutoff && sigma[index_ns] > 0)
    {
        C.setZero();
        sigma.setZero();
        material_state_variables.setTensileStress(true);
    }
}

template class LinearElasticIsotropic<2>;
template class LinearElasticIsotropic<3>;


}   // namespace Fracture
}   // namespace MaterialLib
