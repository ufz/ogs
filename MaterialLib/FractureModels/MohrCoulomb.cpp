/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
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
        Eigen::Ref<Eigen::MatrixXd> Kep)
{
    if (DisplacementDim == 3)
    {
        OGS_FATAL("MohrCoulomb fracture model does not support 3D case.");
        return;
    }
    MaterialPropertyValues const mat(_mp, t, x);
    Eigen::VectorXd const dw = w - w_prev;

    Eigen::MatrixXd Ke(2,2);
    Ke.setZero();
    Ke(0,0) = mat.Ks;
    Ke(1,1) = mat.Kn;

    sigma.noalias() = sigma_prev + Ke * dw;

    // if opening
    if (sigma[1] > 0)
    {
        Kep.setZero();
        sigma.setZero();
        return;
    }

    // check shear yield function (Fs)
    double const Fs = std::abs(sigma[0]) + sigma[1] * tan(mat.phi) - mat.c;
    if (Fs < .0)
    {
        Kep = Ke;
        return;
    }

    Eigen::VectorXd dFs_dS(2);
    dFs_dS[0] = MathLib::sgn(sigma[0]);
    dFs_dS[1] = tan(mat.phi);

    // plastic potential function: Qs = |tau| + Sn * tan da
    Eigen::VectorXd dQs_dS(2);
    dQs_dS[0] = MathLib::sgn(sigma[0]);
    dQs_dS[1] = tan(mat.psi);

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
