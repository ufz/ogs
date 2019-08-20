/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MathLib/KelvinVector.h"
#include "NumLib/NewtonRaphson.h"
#include "ParameterLib/Parameter.h"

#include "MechanicsBase.h"

namespace MaterialLib
{
namespace Solids
{
namespace Lubby2
{
//
// Variables specific to the material model.
//
struct Lubby2MaterialProperties
{
    using P = ParameterLib::Parameter<double>;
    Lubby2MaterialProperties(P const& GK0_,
                             P const& GM0_,
                             P const& KM0_,
                             P const& etaK0_,
                             P const& etaM0_,
                             P const& mK_,
                             P const& mvK_,
                             P const& mvM_)
        : GK0(GK0_),
          GM0(GM0_),
          KM0(KM0_),
          etaK0(etaK0_),
          etaM0(etaM0_),
          mK(mK_),
          mvK(mvK_),
          mvM(mvM_)
    {
    }

    // basic material parameters
    P const& GK0;
    P const& GM0;
    P const& KM0;
    P const& etaK0;
    P const& etaM0;
    P const& mK;
    P const& mvK;
    P const& mvM;
};

namespace detail
{
template <int DisplacementDim>
struct LocalLubby2Properties
{
    LocalLubby2Properties(double const t,
                          ParameterLib::SpatialPosition const& x,
                          Lubby2MaterialProperties const& mp)
        : GM0(mp.GM0(t, x)[0]),
          KM0(mp.KM0(t, x)[0]),
          GK0(mp.GK0(t, x)[0]),
          etaK0(mp.etaK0(t, x)[0]),
          etaM0(mp.etaM0(t, x)[0]),
          mK(mp.mK(t, x)[0]),
          mvK(mp.mvK(t, x)[0]),
          mvM(mp.mvM(t, x)[0])
    {
    }

    void update(double const s_eff)
    {
        double const GM0_s_eff = GM0 * s_eff;
        GK = GK0 * std::exp(mK * GM0_s_eff);
        etaK = etaK0 * std::exp(mvK * GM0_s_eff);
        etaM = etaM0 * std::exp(mvM * GM0_s_eff);
    }

    double const GM0;
    double const KM0;
    double const GK0;
    double const etaK0;
    double const etaM0;
    double const mK;
    double const mvK;
    double const mvM;

    // Solution dependent values.
    double GK = std::numeric_limits<double>::quiet_NaN();
    double etaK = std::numeric_limits<double>::quiet_NaN();
    double etaM = std::numeric_limits<double>::quiet_NaN();
};
}  // namespace detail

template <int DisplacementDim>
class Lubby2 final : public MechanicsBase<DisplacementDim>
{
public:
    struct MaterialStateVariables
        : public MechanicsBase<DisplacementDim>::MaterialStateVariables
    {
        MaterialStateVariables()
        {
            // Previous time step values are not initialized and are set later.
            eps_K_t.resize(kelvin_vector_size);
            eps_M_t.resize(kelvin_vector_size);

            // Initialize current time step values
            eps_K_j.setZero(kelvin_vector_size);
            eps_M_j.setZero(kelvin_vector_size);
        }

        void setInitialConditions()
        {
            eps_K_j = eps_K_t;
            eps_M_j = eps_M_t;
        }

        void pushBackState() override
        {
            eps_K_t = eps_K_j;
            eps_M_t = eps_M_j;
        }

        using KelvinVector =
            MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;
        using KelvinMatrix =
            MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>;
        /// Deviatoric strain in the viscous kelvin element during the current
        /// iteration
        KelvinVector eps_K_t;
        KelvinVector eps_K_j;
        /// Deviatoric strain in the viscous maxwell element during the current
        /// iteration
        KelvinVector eps_M_t;
        KelvinVector eps_M_j;

        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    };

    std::unique_ptr<
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables>
    createMaterialStateVariables() const override
    {
        return std::unique_ptr<
            typename MechanicsBase<DisplacementDim>::MaterialStateVariables>{
            new MaterialStateVariables};
    }

    double getBulkModulus(double const t,
                          ParameterLib::SpatialPosition const& x) const override
    {
        return _mp.KM0(t, x)[0];
    }

public:
    static int const kelvin_vector_size =
        MathLib::KelvinVector::size<DisplacementDim>();
    using KelvinVector =
        MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;
    using KelvinMatrix =
        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>;

    static int const JacobianResidualSize =
        3 * kelvin_vector_size;  // Three is the number of components in the
                                 // jacobian/residual, not the space dimension.
    using ResidualVector = Eigen::Matrix<double, JacobianResidualSize, 1>;
    using JacobianMatrix = Eigen::Matrix<double,
                                         JacobianResidualSize,
                                         JacobianResidualSize,
                                         Eigen::RowMajor>;

public:
    explicit Lubby2(
        NumLib::NewtonRaphsonSolverParameters nonlinear_solver_parameters,
        Lubby2MaterialProperties& material_properties)
        : _nonlinear_solver_parameters(std::move(nonlinear_solver_parameters)),
          _mp(material_properties)
    {
    }

    double computeFreeEnergyDensity(
        double const t,
        ParameterLib::SpatialPosition const& x,
        double const dt,
        KelvinVector const& eps,
        KelvinVector const& sigma,
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables const&
            material_state_variables) const override
    {
        assert(dynamic_cast<MaterialStateVariables const*>(
                   &material_state_variables) != nullptr);
        MaterialStateVariables state(static_cast<MaterialStateVariables const&>(
            material_state_variables));

        auto const& eps_K = state.eps_K_j;
        auto const& eps_K_prev = state.eps_K_t;

        auto const& eps_M = state.eps_M_j;

        auto local_lubby2_properties =
            detail::LocalLubby2Properties<DisplacementDim>{t, x, _mp};

        // calculation of deviatoric parts
        using Invariants =
            MathLib::KelvinVector::Invariants<kelvin_vector_size>;
        auto const& P_dev = Invariants::deviatoric_projection;
        KelvinVector const epsd_i = P_dev * eps;

        // initial guess as elastic predictor
        KelvinVector sigd_j = 2.0 * (epsd_i - state.eps_M_t - state.eps_K_t);

        // Calculate effective stress and update material properties
        double sig_eff = Invariants::equivalentStress(sigd_j);
        local_lubby2_properties.update(sig_eff);

        auto const& eta_K = local_lubby2_properties.etaK;

        // This is specific to the backward Euler time scheme and needs to be
        // updated if other time schemes are used.
        return (eps - eps_K - eps_M).dot(sigma) / 2 +
               eps_K.dot(sigma - eta_K * (eps_K - eps_K_prev) / dt) / 2;
    }

    boost::optional<std::tuple<KelvinVector,
                               std::unique_ptr<typename MechanicsBase<
                                   DisplacementDim>::MaterialStateVariables>,
                               KelvinMatrix>>
    integrateStress(
        double const t, ParameterLib::SpatialPosition const& x, double const dt,
        KelvinVector const& eps_prev, KelvinVector const& eps,
        KelvinVector const& sigma_prev,
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables const&
            material_state_variables,
        double const T) const override;

private:
    /// Calculates the 18x1 residual vector.
    void calculateResidualBurgers(
        double const dt,
        const KelvinVector& strain_curr,
        const KelvinVector& strain_t,
        const KelvinVector& stress_curr,
        const KelvinVector& stress_t,
        KelvinVector& strain_Kel_curr,
        const KelvinVector& strain_Kel_t,
        KelvinVector& strain_Max_curr,
        const KelvinVector& strain_Max_t,
        ResidualVector& res,
        detail::LocalLubby2Properties<DisplacementDim> const& properties) const;

    /// Calculates the 18x18 Jacobian.
    void calculateJacobianBurgers(
        double const t,
        ParameterLib::SpatialPosition const& x,
        double const dt,
        JacobianMatrix& Jac,
        double s_eff,
        const KelvinVector& sig_i,
        const KelvinVector& eps_K_i,
        detail::LocalLubby2Properties<DisplacementDim> const& properties) const;

private:
    NumLib::NewtonRaphsonSolverParameters const _nonlinear_solver_parameters;
    Lubby2MaterialProperties _mp;
};

extern template class Lubby2<2>;
extern template class Lubby2<3>;

}  // namespace Lubby2
}  // namespace Solids
}  // namespace MaterialLib
