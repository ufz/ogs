/**
 * \file
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * Implementation of Ehler's single-surface model.
 * see Ehler's paper "A single-surface yield function for geomaterials" for more
 * details. \cite Ehlers1995
 *
 * Refer to "Single-surface benchmark of OpenGeoSys documentation
 * (https://docs.opengeosys.org/docs/benchmarks/small-deformations/mechanics-plasticity-single-surface)"
 * for more details for the tests.
 */

#pragma once

#ifndef NDEBUG
#include <ostream>
#endif

#include "BaseLib/Error.h"
#include "MathLib/KelvinVector.h"
#include "NumLib/NewtonRaphson.h"
#include "ProcessLib/Parameter/Parameter.h"

#include "MechanicsBase.h"

namespace MaterialLib
{
namespace Solids
{
namespace Ehlers
{
/// material parameters in relation to Ehler's single-surface model see Ehler's
/// paper "A single-surface yield function for geomaterials" for more details
/// \cite Ehlers1995.
struct MaterialPropertiesParameters
{
    using P = ProcessLib::Parameter<double>;

    MaterialPropertiesParameters(P const& G_, P const& K_, P const& alpha_,
                                 P const& beta_, P const& gamma_,
                                 P const& delta_, P const& epsilon_,
                                 P const& m_, P const& alpha_p_,
                                 P const& beta_p_, P const& gamma_p_,
                                 P const& delta_p_, P const& epsilon_p_,
                                 P const& m_p_, P const& kappa_,
                                 P const& hardening_coefficient_)
        : G(G_),
          K(K_),
          alpha(alpha_),
          beta(beta_),
          gamma(gamma_),
          delta(delta_),
          epsilon(epsilon_),
          m(m_),
          alpha_p(alpha_p_),
          beta_p(beta_p_),
          gamma_p(gamma_p_),
          delta_p(delta_p_),
          epsilon_p(epsilon_p_),
          m_p(m_p_),
          kappa(kappa_),
          hardening_coefficient(hardening_coefficient_)
    {
    }
    // basic material parameters
    P const& G;  ///< shear modulus
    P const& K;  ///< bulk modulus

    P const& alpha;    ///< material dependent parameter in relation to Ehlers
                       ///< model, refer to \cite Ehlers1995 .
    P const& beta;     ///< \copydoc alpha
    P const& gamma;    ///< \copydoc alpha
    P const& delta;    ///< \copydoc alpha
    P const& epsilon;  ///< \copydoc alpha
    P const& m;        ///< \copydoc alpha

    P const& alpha_p;    ///< \copydoc alpha
    P const& beta_p;     ///< \copydoc alpha
    P const& gamma_p;    ///< \copydoc alpha
    P const& delta_p;    ///< \copydoc alpha
    P const& epsilon_p;  ///< \copydoc alpha
    P const& m_p;        ///< \copydoc alpha

    P const& kappa;  ///< hardening parameter
    P const& hardening_coefficient;
};

struct DamagePropertiesParameters
{
    using P = ProcessLib::Parameter<double>;
    P const& alpha_d;
    P const& beta_d;
    P const& h_d;
};

/// Evaluated MaterialPropertiesParameters container, see its documentation for
/// details.
struct MaterialProperties final
{
    MaterialProperties(double const t, ProcessLib::SpatialPosition const& x,
                       MaterialPropertiesParameters const& mp)
        : G(mp.G(t, x)[0]),
          K(mp.K(t, x)[0]),
          alpha(mp.alpha(t, x)[0]),
          beta(mp.beta(t, x)[0]),
          gamma(mp.gamma(t, x)[0]),
          delta(mp.delta(t, x)[0]),
          epsilon(mp.epsilon(t, x)[0]),
          m(mp.m(t, x)[0]),
          alpha_p(mp.alpha_p(t, x)[0]),
          beta_p(mp.beta_p(t, x)[0]),
          gamma_p(mp.gamma_p(t, x)[0]),
          delta_p(mp.delta_p(t, x)[0]),
          epsilon_p(mp.epsilon_p(t, x)[0]),
          m_p(mp.m_p(t, x)[0]),
          kappa(mp.kappa(t, x)[0]),
          hardening_coefficient(mp.hardening_coefficient(t, x)[0])
    {
    }
    double const G;
    double const K;

    double const alpha;
    double const beta;
    double const gamma;
    double const delta;
    double const epsilon;
    double const m;

    double const alpha_p;
    double const beta_p;
    double const gamma_p;
    double const delta_p;
    double const epsilon_p;
    double const m_p;

    double const kappa;
    double const hardening_coefficient;
};

/// Evaluated DamagePropertiesParameters container, see its documentation for
/// details.
struct DamageProperties
{
    DamageProperties(double const t,
                     ProcessLib::SpatialPosition const& x,
                     DamagePropertiesParameters const& dp)
        : alpha_d(dp.alpha_d(t, x)[0]),
          beta_d(dp.beta_d(t, x)[0]),
          h_d(dp.h_d(t, x)[0])
    {
    }
    double const alpha_d;
    double const beta_d;
    double const h_d;
};

template <typename KelvinVector>
struct PlasticStrain final
{
    PlasticStrain() : D(KelvinVector::Zero()) {}
    PlasticStrain(KelvinVector eps_p_D_, double const eps_p_V_,
                  double const eps_p_eff_)
        : D(std::move(eps_p_D_)), V(eps_p_V_), eff(eps_p_eff_){};

    KelvinVector D;  ///< deviatoric plastic strain
    double V = 0;    ///< volumetric strain
    double eff = 0;  ///< effective plastic strain

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

class Damage final
{
public:
    Damage() = default;
    Damage(double const kappa_d, double const value)
        : _kappa_d(kappa_d), _value(value)
    {
    }

    double kappa_d() const { return _kappa_d; }
    double value() const { return _value; }

private:
    double _kappa_d = 0;  ///< damage driving variable
    double _value = 0;    ///< isotropic damage variable
};

template <int DisplacementDim>
struct StateVariables
    : public MechanicsBase<DisplacementDim>::MaterialStateVariables
{
    StateVariables& operator=(StateVariables const&) = default;
    typename MechanicsBase<DisplacementDim>::MaterialStateVariables& operator=(
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables const&
            state) noexcept override
    {
        assert(dynamic_cast<StateVariables const*>(&state) != nullptr);
        return operator=(static_cast<StateVariables const&>(state));
    }

    void setInitialConditions()
    {
        eps_p = eps_p_prev;
        damage = damage_prev;
    }

    void pushBackState() override
    {
        eps_p_prev = eps_p;
        damage_prev = damage;
    }

    using KelvinVector =
        MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;

    PlasticStrain<KelvinVector> eps_p;  ///< plastic part of the state.
    Damage damage;                      ///< damage part of the state.

    // Initial values from previous timestep
    PlasticStrain<KelvinVector> eps_p_prev;  ///< \copydoc eps_p
    Damage damage_prev;                      ///< \copydoc damage

#ifndef NDEBUG
    friend std::ostream& operator<<(
        std::ostream& os, StateVariables<DisplacementDim> const& m)
    {
        os << "State:\n"
           << "eps_p_D: " << m.eps_p.D << "\n"
           << "eps_p_eff: " << m.eps_p.eff << "\n"
           << "kappa_d: " << m.damage.kappa_d() << "\n"
           << "damage: " << m.damage.value() << "\n"
           << "eps_p_D_prev: " << m.eps_p_prev.D << "\n"
           << "eps_p_eff_prev: " << m.eps_p_prev.eff << "\n"
           << "kappa_d_prev: " << m.damage_prev.kappa_d() << "\n"
           << "damage_prev: " << m.damage_prev.value() << "\n"
           << "lambda: " << m.lambda << "\n";
        return os;
    }
#endif  // NDEBUG

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

template <int DisplacementDim>
class SolidEhlers final : public MechanicsBase<DisplacementDim>
{
public:
    static int const KelvinVectorSize =
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value;
    static int const JacobianResidualSize =
        2 * KelvinVectorSize + 3;  // 2 is the number of components in the
                                   // jacobian/residual, not the space
                                   // dimension. And 3 is for additional
                                   // variables.
    using ResidualVectorType = Eigen::Matrix<double, JacobianResidualSize, 1>;
    using JacobianMatrix = Eigen::Matrix<double, JacobianResidualSize,
                                         JacobianResidualSize, Eigen::RowMajor>;

public:
    std::unique_ptr<
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables>
    createMaterialStateVariables() const override
    {
        return std::make_unique<StateVariables<DisplacementDim>>();
    }

public:
    using KelvinVector =
        MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;
    using KelvinMatrix =
        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>;

public:
    explicit SolidEhlers(
        NumLib::NewtonRaphsonSolverParameters nonlinear_solver_parameters,
        MaterialPropertiesParameters material_properties,
        std::unique_ptr<DamagePropertiesParameters>&& damage_properties)
        : _nonlinear_solver_parameters(std::move(nonlinear_solver_parameters)),
          _mp(std::move(material_properties)),
          _damage_properties(std::move(damage_properties))
    {
    }

    double computeFreeEnergyDensity(
        double const /*t*/,
        ProcessLib::SpatialPosition const& /*x*/,
        double const /*dt*/,
        KelvinVector const& eps,
        KelvinVector const& sigma,
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables const&
            material_state_variables) const override
    {
        assert(dynamic_cast<StateVariables<DisplacementDim> const*>(
                   &material_state_variables) != nullptr);

        auto const& eps_p = static_cast<StateVariables<DisplacementDim> const&>(
                                material_state_variables)
                                .eps_p;
        using Invariants = MathLib::KelvinVector::Invariants<KelvinVectorSize>;
        auto const& identity2 = Invariants::identity2;
        return (eps - eps_p.D - eps_p.V / 3 * identity2).dot(sigma) / 2;
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
        double const T) const override;

    std::vector<typename MechanicsBase<DisplacementDim>::InternalVariable>
    getInternalVariables() const override;

    MaterialProperties evaluatedMaterialProperties(
        double const t, ProcessLib::SpatialPosition const& x) const
    {
        return MaterialProperties(t, x, _mp);
    }

    DamageProperties evaluatedDamageProperties(
        double const t, ProcessLib::SpatialPosition const& x) const
    {
        return DamageProperties(t, x, *_damage_properties);
    }

private:
    NumLib::NewtonRaphsonSolverParameters const _nonlinear_solver_parameters;

    MaterialPropertiesParameters _mp;
    std::unique_ptr<DamagePropertiesParameters> _damage_properties;
};

extern template class SolidEhlers<2>;
extern template class SolidEhlers<3>;
}  // namespace Ehlers
}  // namespace Solids
}  // namespace MaterialLib
