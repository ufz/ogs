/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MohrCoulomb.h"

#include "BaseLib/Error.h"
#include "MathLib/MathTools.h"

namespace MaterialLib
{
namespace Fracture
{

namespace
{

struct MaterialPropertyValues
{
    double Kn = 0.0;
    double Ks = 0.0;
    double phi = 0.0; // friction angle
    double psi = 0.0; // dilation angle
    double c = 0.0;

    template <typename MaterialProperties>
    MaterialPropertyValues(
            MaterialProperties const& mp,
            double const t,
            ProcessLib::SpatialPosition const& x)
    {
        Kn = mp.normal_stiffness(t,x)[0];
        Ks = mp.shear_stiffness(t,x)[0];
        phi = MathLib::to_radians(mp.friction_angle(t,x)[0]);
        psi = MathLib::to_radians(mp.dilatancy_angle(t,x)[0]);
        c = mp.cohesion(t,x)[0];
    }
};

} // no namespace

template <int DisplacementDim>
void MohrCoulomb<DisplacementDim>::computeConstitutiveRelation(
        double const t,
        ProcessLib::SpatialPosition const& x,
        Eigen::Ref<Eigen::VectorXd const> w_prev,
        Eigen::Ref<Eigen::VectorXd const> w,
        Eigen::Ref<Eigen::VectorXd const> sigma_prev,
        Eigen::Ref<Eigen::VectorXd> sigma,
        Eigen::Ref<Eigen::MatrixXd> Kep,
        typename FractureModelBase<DisplacementDim>::MaterialStateVariables&
        material_state_variables)
{
    material_state_variables.reset();

    MaterialPropertyValues const mat(_mp, t, x);
    Eigen::VectorXd const dw = w - w_prev;

    const int index_ns = DisplacementDim - 1;
    Eigen::MatrixXd Ke = Eigen::MatrixXd::Zero(DisplacementDim,DisplacementDim);
    for (int i=0; i<index_ns; i++)
        Ke(i,i) = mat.Ks;
    Ke(index_ns, index_ns) = mat.Kn;

    sigma.noalias() = sigma_prev + Ke * dw;
    double const tau = (DisplacementDim==2) ? sigma[0] : std::sqrt(sigma[0]*sigma[0] + sigma[1]*sigma[1]);
    double const sigma_n = sigma[index_ns];

    // if opening
    if (sigma_n > 0)
    {
        Kep.setZero();
        sigma.setZero();
        material_state_variables.setTensileStress(true);
        return;
    }

    // check shear yield function (Fs)
    double const Fs = std::abs(tau) + sigma_n * std::tan(mat.phi) - mat.c;

    material_state_variables.setShearYieldFunctionValue(Fs);
    if (Fs < .0)
    {
        Kep = Ke;
        return;
    }

    Eigen::VectorXd dFs_dS(DisplacementDim);
    // plastic potential function: Qs = |tau| + Sn * tan da
    Eigen::VectorXd dQs_dS(DisplacementDim);

    if (DisplacementDim == 2)
    {
        dFs_dS[0] = boost::math::sign(tau);
        dQs_dS[0] = boost::math::sign(tau);
    }
    else
    {
        for (int i=0; i<index_ns; i++)
            dFs_dS[i] = boost::math::sign(tau) / tau * sigma[i];

        for (int i=0; i<index_ns; i++)
            dQs_dS[i] = boost::math::sign(tau) / tau * sigma[i];
    }
    dFs_dS[index_ns] = std::tan(mat.phi);
    dQs_dS[index_ns] = std::tan(mat.psi);

    // plastic multiplier
    Eigen::RowVectorXd const A = dFs_dS.transpose() * Ke / (dFs_dS.transpose() * Ke * dQs_dS);
    double const d_eta = A * dw;

    // plastic part of the dispalcement
    Eigen::VectorXd const dwp = dQs_dS * d_eta;

    // correct stress
    sigma.noalias() = sigma_prev + Ke * (dw - dwp);

    // Kep
    Kep = Ke - Ke * dQs_dS * A;
}

template class MohrCoulomb<2>;
template class MohrCoulomb<3>;

}   // namespace Fracture
}  // namespace MaterialLib
