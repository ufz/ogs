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

#ifndef MATERIALLIB_SOLIDMODELS_EHLERS_H_
#define MATERIALLIB_SOLIDMODELS_EHLERS_H_

#include <cfloat>
#include <memory>
#ifndef NDEBUG
#include <ostream>
#endif

#include <Eigen/Dense>
#include <logog/include/logog.hpp>

#include "BaseLib/Error.h"
#include "KelvinVector.h"
#include "MechanicsBase.h"

namespace MaterialLib
{
namespace Solids
{

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
    //
    // Variables specific to the material model.
    //
    struct MaterialProperties
    {
        using P = ProcessLib::Parameter<double>;

        using KelvinVector = ProcessLib::KelvinVectorType<DisplacementDim>;
        using KelvinMatrix = ProcessLib::KelvinMatrixType<DisplacementDim>;

        /// material parameters in relation to Ehler's single-surface model
        /// see Ehler's paper "A single-surface yield function for geomaterials"
        /// for more details.
        MaterialProperties(P const& G_, P const& K_, P const& alpha_,
                           P const& beta_, P const& gamma_, P const& delta_,
                           P const& epsilon_, P const& m_, P const& alpha_p_,
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
        // Drucker-Prager: Import kappa and beta in terms of Drucker-Prager
        // criterion solution dependent values
        double k;

        void calculateIsotropicHardening(double const t,
                                         ProcessLib::SpatialPosition const& x,
                                         const double e_pv_curr);
    };

    struct MaterialStateVariables
        : public MechanicsBase<DisplacementDim>::MaterialStateVariables
    {
        MaterialStateVariables()
            : eps_p_D(KelvinVector::Zero()), eps_p_D_prev(KelvinVector::Zero())
        {
        }

        void loadState()
        {
            eps_p_D = eps_p_D_prev;
            eps_p_V = eps_p_V_prev;
            eps_p_eff = eps_p_eff_prev;
            lambda = 0;
        }

        void pushBackState() override
        {
            eps_p_D_prev = eps_p_D;
            eps_p_V_prev = eps_p_V;
            eps_p_eff_prev = eps_p_eff;  // effective part of trace(eps_p)
            lambda = 0;
        }

        using KelvinVector = ProcessLib::KelvinVectorType<DisplacementDim>;

        KelvinVector eps_p_D;  ///< deviatoric plastic strain
        double eps_p_V = 0;    ///< volumetric strain
        double eps_p_eff = 0;  ///< effective plastic strain

        // Initial values from previous timestep
        KelvinVector eps_p_D_prev;  ///< \copydoc eps_p_D
        double eps_p_V_prev = 0;    ///< \copydoc eps_p_V
        double eps_p_eff_prev = 0;  ///< \copydoc eps_p_eff

        double lambda = 0;  ///< plastic multiplier

#ifndef NDEBUG
        friend std::ostream& operator<<(std::ostream& os,
                                        MaterialStateVariables const& m)
        {
            os << "State:\n"
               << "eps_p_D: " << m.eps_p_D << "\n"
               << "eps_p_eff: " << m.eps_p_eff << "\n"
               << "eps_p_D_prev: " << m.eps_p_D_prev << "\n"
               << "eps_p_eff_prev: " << m.eps_p_eff_prev << "\n"
               << "lambda: " << m.lambda << "\n";
            return os;
        }
#endif  // NDEBUG
    };

    std::unique_ptr<
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables>
    createMaterialStateVariables() override
    {
        return std::unique_ptr<
            typename MechanicsBase<DisplacementDim>::MaterialStateVariables>{
            new MaterialStateVariables};
    }

public:
    using KelvinVector = ProcessLib::KelvinVectorType<DisplacementDim>;
    using KelvinMatrix = ProcessLib::KelvinMatrixType<DisplacementDim>;

public:
    explicit SolidEhlers(MaterialProperties const& material_properties)
        : _mp(material_properties)
    {
    }

    bool computeConstitutiveRelation(
        double const t,
        ProcessLib::SpatialPosition const& x,
        double const dt,
        KelvinVector const& eps_prev,
        KelvinVector const& eps,
        KelvinVector const& sigma_prev,
        KelvinVector& sigma,
        KelvinMatrix& C,
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables&
            material_state_variables) override;

private:
    MaterialProperties _mp;
};

}  // namespace Solids
}  // namespace MaterialLib
#include "Ehlers-impl.h"

#endif  // MATERIALLIB_SOLIDMODELS_EHLERS_H_
