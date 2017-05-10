/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/*
 * Implementation of Ehler's single-surface model.
 * see Ehler's paper "A single-surface yield function for geomaterials" for more
 * details. \cite{Ehlers1995}
 *
 * Refer to "Single-surface benchmark of OpenGeoSys documentation
 * (https://docs.opengeosys.org/docs/benchmarks/small-deformations/mechanics-plasticity-single-surface)"
 * for more details for the tests.
 */

#pragma once

#include <cfloat>
#include <memory>
#ifndef NDEBUG
#include <ostream>
#endif

#include <Eigen/Dense>
#include <logog/include/logog.hpp>
#include <utility>

#include "BaseLib/Error.h"
#include "NumLib/NewtonRaphson.h"

#include "KelvinVector.h"
#include "MechanicsBase.h"

namespace MaterialLib
{
namespace Solids
{
namespace Ehlers
{
//
// Variables specific to the material model.
//
struct MaterialPropertiesParameters
{
    using P = ProcessLib::Parameter<double>;

    /// material parameters in relation to Ehler's single-surface model
    /// see Ehler's paper "A single-surface yield function for geomaterials"
    /// for more details.
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

    P const& alpha;
    P const& beta;
    P const& gamma;
    P const& delta;
    P const& epsilon;
    P const& m;

    P const& alpha_p;
    P const& beta_p;
    P const& gamma_p;
    P const& delta_p;
    P const& epsilon_p;
    P const& m_p;

    P const& kappa;
    P const& hardening_coefficient;
};

struct DamagePropertiesParameters
{
    using P = ProcessLib::Parameter<double>;
    P const& alpha_d;
    P const& beta_d;
    P const& h_d;
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
};

struct Damage final
{
    Damage() = default;
    Damage(double const kappa_d_, double const damage_)
        : kappa_d(kappa_d_), damage(damage_){};
    double kappa_d = 0;  ///< damage driving variable
    double damage = 0;   ///< isotropic damage variable
};

template <int DisplacementDim>
struct MaterialStateVariables
    : public MechanicsBase<DisplacementDim>::MaterialStateVariables
{
    MaterialStateVariables& operator=(MaterialStateVariables const&) = default;
    typename MechanicsBase<DisplacementDim>::MaterialStateVariables& operator=(
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables const&
            state) noexcept override
    {
        assert(dynamic_cast<MaterialStateVariables const*>(&state) != nullptr);
        return operator=(static_cast<MaterialStateVariables const&>(state));
    }

    void setInitialConditions()
    {
        eps_p = eps_p_prev;
        lambda = 0;
        damage = damage_prev;
    }

    void pushBackState() override
    {
        eps_p_prev = eps_p;
        lambda = 0;
        damage_prev = damage;
    }

    using KelvinVector = ProcessLib::KelvinVectorType<DisplacementDim>;

    PlasticStrain<KelvinVector> eps_p;  ///< plastic part of the state.
    Damage damage;                      ///< damage part of the state.

    // Initial values from previous timestep
    PlasticStrain<KelvinVector> eps_p_prev;  ///< \copydoc eps_p
    double lambda = 0;                       ///< plastic multiplier
    Damage damage_prev;                      ///< \copydoc damage

#ifndef NDEBUG
    friend std::ostream& operator<<(
        std::ostream& os, MaterialStateVariables<DisplacementDim> const& m)
    {
        os << "State:\n"
           << "eps_p_D: " << m.eps_p.D << "\n"
           << "eps_p_eff: " << m.eps_p.eff << "\n"
           << "kappa_d: " << m.damage.kappa_d << "\n"
           << "damage: " << m.damage.damage << "\n"
           << "eps_p_D_prev: " << m.eps_p_prev.D << "\n"
           << "eps_p_eff_prev: " << m.eps_p_prev.eff << "\n"
           << "kappa_d_prev: " << m.damage_prev.kappa_d << "\n"
           << "damage_prev: " << m.damage_prev.damage << "\n"
           << "lambda: " << m.lambda << "\n";
        return os;
    }
#endif  // NDEBUG
};

template <int DisplacementDim>
class SolidEhlers final : public MechanicsBase<DisplacementDim>
{
public:
    static int const KelvinVectorSize =
        ProcessLib::KelvinVectorDimensions<DisplacementDim>::value;
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
    createMaterialStateVariables() override
    {
        return std::unique_ptr<
            typename MechanicsBase<DisplacementDim>::MaterialStateVariables>{
            new MaterialStateVariables<DisplacementDim>};
    }

public:
    using KelvinVector = ProcessLib::KelvinVectorType<DisplacementDim>;
    using KelvinMatrix = ProcessLib::KelvinMatrixType<DisplacementDim>;

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

    std::tuple<KelvinVector, std::unique_ptr<typename MechanicsBase<
                                 DisplacementDim>::MaterialStateVariables>,
               KelvinMatrix>
    integrateStress(
        double const t,
        ProcessLib::SpatialPosition const& x,
        double const dt,
        KelvinVector const& eps_prev,
        KelvinVector const& eps,
        KelvinVector const& sigma_prev,
        KelvinVector const& sigma,
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables const&
            material_state_variables) override;

private:
    NumLib::NewtonRaphsonSolverParameters const _nonlinear_solver_parameters;

    MaterialPropertiesParameters _mp;
    std::unique_ptr<DamagePropertiesParameters> _damage_properties;
};

}  // namespace Ehlers
}  // namespace Solids
}  // namespace MaterialLib
#include "Ehlers-impl.h"
