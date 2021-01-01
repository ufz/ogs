/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "Coulomb.h"

#include <boost/math/constants/constants.hpp>

#include "BaseLib/Error.h"
#include "LogPenalty.h"
#include "MathLib/MathTools.h"
#include "NumLib/Exceptions.h"

namespace MaterialLib
{
namespace Fracture
{
namespace Coulomb
{

struct MaterialPropertyValues
{
    double Kn = 0.0;
    double Ks = 0.0;
    double phi = 0.0; // friction angle
    double psi = 0.0; // dilation angle
    double c = 0.0;

    template <typename MaterialProperties>
    MaterialPropertyValues(MaterialProperties const& mp,
                           double const t,
                           ParameterLib::SpatialPosition const& x)
    {
        Kn = mp.normal_stiffness(t,x)[0];
        Ks = mp.shear_stiffness(t,x)[0];
        auto constexpr degree =
            boost::math::constants::degree<double>();  // pi/180
        phi = mp.friction_angle(t, x)[0] * degree;
        psi = mp.dilatancy_angle(t, x)[0] * degree;
        c = mp.cohesion(t,x)[0];
    }
};

template <int DisplacementDim>
void Coulomb<DisplacementDim>::computeConstitutiveRelation(
    double const t,
    ParameterLib::SpatialPosition const& x,
    double const aperture0,
    Eigen::Ref<Eigen::VectorXd const>
    /*sigma0*/,
    Eigen::Ref<Eigen::VectorXd const>
        w_prev,
    Eigen::Ref<Eigen::VectorXd const>
        w,
    Eigen::Ref<Eigen::VectorXd const>
        sigma_prev,
    Eigen::Ref<Eigen::VectorXd>
        sigma,
    Eigen::Ref<Eigen::MatrixXd>
        Kep,
    typename FractureModelBase<DisplacementDim>::MaterialStateVariables&
        material_state_variables)
{
    assert(dynamic_cast<StateVariables<DisplacementDim> const*>(
               &material_state_variables) != nullptr);

    StateVariables<DisplacementDim>& state =
        static_cast<StateVariables<DisplacementDim>&>(material_state_variables);

    MaterialPropertyValues const mat(_mp, t, x);

    const int index_ns = DisplacementDim - 1;
    double const aperture = w[index_ns] + aperture0;

    Eigen::MatrixXd Ke;
    {  // Elastic tangent stiffness
        Ke = Eigen::MatrixXd::Zero(DisplacementDim, DisplacementDim);
        for (int i = 0; i < index_ns; i++)
        {
            Ke(i, i) = mat.Ks;
        }

        Ke(index_ns, index_ns) = mat.Kn;
    }

    // Total plastic aperture compression
    // NOTE: Initial condition sigma0 seems to be associated with an initial
    // condition of the w0 = 0. Therefore the initial state is not associated
    // with a plastic aperture change.
    {  // Exact elastic predictor
        sigma.noalias() = Ke * (w - w_prev);

        sigma.coeffRef(index_ns) *=
            logPenaltyDerivative(aperture0, aperture, _penalty_aperture_cutoff);
        sigma.noalias() += sigma_prev;
    }

    // correction for an opening fracture
    if (_tension_cutoff && sigma[DisplacementDim - 1] >= 0)
    {
        Kep.setZero();
        sigma.setZero();
        state.w_p = w;
        material_state_variables.setTensileStress(true);
        return;
    }

    auto yield_function = [&mat](Eigen::VectorXd const& s) {
        double const sigma_n = s[DisplacementDim - 1];
        Eigen::VectorXd const sigma_s = s.head(DisplacementDim - 1);
        double const mag_tau = sigma_s.norm();  // magnitude
        return mag_tau + sigma_n * std::tan(mat.phi) - mat.c;
    };

    {  // Exit if still in elastic range by checking the shear yield function.
        double const Fs = yield_function(sigma);
        material_state_variables.setShearYieldFunctionValue(Fs);
        if (Fs < .0)
        {
            Kep = Ke;
            Kep(index_ns, index_ns) *= logPenaltyDerivative(
                aperture0, aperture, _penalty_aperture_cutoff);
            return;
        }
    }

    auto yield_function_derivative = [&mat](Eigen::VectorXd const& s) {
        Eigen::Matrix<double, DisplacementDim, 1> dFs_dS;
        dFs_dS.template head<DisplacementDim - 1>().noalias() =
            s.template head<DisplacementDim - 1>().normalized();
        dFs_dS.coeffRef(DisplacementDim - 1) = std::tan(mat.phi);
        return dFs_dS;
    };

    // plastic potential function: Qs = |tau| + Sn * tan da
    auto plastic_potential_derivative = [&mat](Eigen::VectorXd const& s) {
        Eigen::Matrix<double, DisplacementDim, 1> dQs_dS;
        dQs_dS.template head<DisplacementDim - 1>().noalias() =
            s.template head<DisplacementDim - 1>().normalized();
        dQs_dS.coeffRef(DisplacementDim - 1) = std::tan(mat.psi);
        return dQs_dS;
    };

    {  // Newton

        Eigen::FullPivLU<Eigen::Matrix<double, 1, 1, Eigen::RowMajor>>
            linear_solver;
        using ResidualVectorType = Eigen::Matrix<double, 1, 1, Eigen::RowMajor>;
        using JacobianMatrix = Eigen::Matrix<double, 1, 1, Eigen::RowMajor>;

        JacobianMatrix jacobian;
        ResidualVectorType solution;
        solution << 0;

        auto const update_residual = [&](ResidualVectorType& residual) {
            residual[0] = yield_function(sigma);
        };

        auto const update_jacobian = [&](JacobianMatrix& jacobian) {
            jacobian(0, 0) = -yield_function_derivative(sigma).transpose() *
                             Ke * plastic_potential_derivative(sigma);
        };

        auto const update_solution = [&](ResidualVectorType const& increment) {
            solution += increment;
            /*DBUG("analytical = {:g}",
                 Fs / (mat.Ks + mat.Kn * std::tan(mat.psi) * std::tan(mat.phi)))
                 */
            state.w_p = state.w_p_prev +
                        solution[0] * plastic_potential_derivative(sigma);

            sigma.noalias() = Ke * (w - w_prev - state.w_p + state.w_p_prev);

            sigma.coeffRef(index_ns) *= logPenaltyDerivative(
                aperture0, aperture, _penalty_aperture_cutoff);
            sigma.noalias() += sigma_prev;
        };

        auto newton_solver =
            NumLib::NewtonRaphson<decltype(linear_solver), JacobianMatrix,
                                  decltype(update_jacobian), ResidualVectorType,
                                  decltype(update_residual),
                                  decltype(update_solution)>(
                linear_solver, update_jacobian, update_residual,
                update_solution, _nonlinear_solver_parameters);

        auto const success_iterations = newton_solver.solve(jacobian);

        if (!success_iterations)
        {
            throw NumLib::AssemblyException(
                "FractureModel/Coulomb local nonlinear solver didn't "
                "converge.");
        }

        // Solution containing lambda is not needed; w_p and sigma already
        // up to date.
    }

    {  // Update material state shear yield function value.
        double const Fs = yield_function(sigma);
        material_state_variables.setShearYieldFunctionValue(Fs);
    }

    Ke(index_ns, index_ns) *=
        logPenaltyDerivative(aperture0, aperture, _penalty_aperture_cutoff);
    Eigen::RowVectorXd const A = yield_function_derivative(sigma).transpose() *
                                 Ke /
                                 (yield_function_derivative(sigma).transpose() *
                                  Ke * plastic_potential_derivative(sigma));
    Kep = Ke - Ke * plastic_potential_derivative(sigma) * A;
}

template class Coulomb<2>;
template class Coulomb<3>;

}  // namespace Coulomb
}  // namespace Fracture
}  // namespace MaterialLib
