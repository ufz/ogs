/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LinearElasticIsotropic.h"

namespace MaterialLib
{
namespace Fracture
{

template <int DisplacementDim>
void LinearElasticIsotropic<DisplacementDim>::computeConstitutiveRelation(
        double const t,
        ProcessLib::SpatialPosition const& x,
        Eigen::Ref<Eigen::VectorXd const> w_prev,
        Eigen::Ref<Eigen::VectorXd const> w,
        Eigen::Ref<Eigen::VectorXd const> sigma_prev,
        Eigen::Ref<Eigen::VectorXd> sigma,
        Eigen::Ref<Eigen::MatrixXd> C,
        typename FractureModelBase<DisplacementDim>::MaterialStateVariables&
        material_state_variables)
{
    material_state_variables.reset();

    const int index_ns = DisplacementDim - 1;
    C.setZero();
    for (int i=0; i<index_ns; i++)
        C(i,i) = _mp.shear_stiffness(t, x)[0];
    C(index_ns, index_ns) = _mp.normal_stiffness(t, x)[0];

    sigma.noalias() = sigma_prev + C * (w - w_prev);

    // correct if fracture is opening
    if (sigma[index_ns] > 0)
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
