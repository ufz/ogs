/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file
 *  Created on July 6, 2018, 9:53 AM
 */

#pragma once

#include <cmath>
#include <memory>
#include <tuple>

#include "LinearElasticIsotropic.h"
#include "MathLib/KelvinVector.h"
#include "NumLib/NewtonRaphson.h"
#include "ParameterLib/Parameter.h"

namespace MaterialLib
{
namespace Solids
{
namespace Creep
{
/**
 * \brief A class for computing the BGRa creep model, which reads
 * \f[
 * \dot {\mathbf{\epsilon}}^{cr}=\sqrt{\frac{3}{2}}A \mathrm{e}^{-\frac{Q}{RT}}
 * \left(\frac{\sigma_{eff}}{\sigma_0}\right)^n\frac{\mathbf{s}}{||\mathbf{s}||}
 * \f]
 * where \f$\sigma_{eff}=\sqrt{\frac{3}{2}}||\mathbf{s}||\f$, \f$A, \sigma_0, n,
 * Q\f$ are parameter, and \f$R\f$ is the gas constant.
 */
template <int DisplacementDim>
class CreepBGRa final : public LinearElasticIsotropic<DisplacementDim>
{
public:
    using LinearElasticIsotropic<DisplacementDim>::KelvinVectorSize;
    using ResidualVectorType = Eigen::Matrix<double, KelvinVectorSize, 1>;
    using JacobianMatrix = Eigen::Matrix<double, KelvinVectorSize,
                                         KelvinVectorSize, Eigen::RowMajor>;

    using KelvinVector =
        MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;
    using KelvinMatrix =
        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>;

    using Parameter = ParameterLib::Parameter<double>;

    using LinearElasticIsotropic<DisplacementDim>::getBulkModulus;

    std::unique_ptr<
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables>
    createMaterialStateVariables() const override
    {
        return LinearElasticIsotropic<
            DisplacementDim>::createMaterialStateVariables();
    }

    CreepBGRa(
        typename LinearElasticIsotropic<DisplacementDim>::MaterialProperties mp,
        NumLib::NewtonRaphsonSolverParameters nonlinear_solver_parameters,
        Parameter const& A, Parameter const& n, Parameter const& sigma_f,
        Parameter const& Q)
        : LinearElasticIsotropic<DisplacementDim>(std::move(mp)),
          _nonlinear_solver_parameters(std::move(nonlinear_solver_parameters)),
          _a(A),
          _n(n),
          _sigma_f(sigma_f),
          _q(Q)
    {
    }

    std::optional<std::tuple<KelvinVector,
                             std::unique_ptr<typename MechanicsBase<
                                 DisplacementDim>::MaterialStateVariables>,
                             KelvinMatrix>>
    integrateStress(
        MaterialPropertyLib::VariableArray const& variable_array_prev,
        MaterialPropertyLib::VariableArray const& variable_array,
        double const t, ParameterLib::SpatialPosition const& x, double const dt,
        KelvinVector const& eps_prev, KelvinVector const& eps,
        KelvinVector const& sigma_prev,
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables const&
            material_state_variables,
        double const T) const override;

    ConstitutiveModel getConstitutiveModel() const override
    {
        return ConstitutiveModel::CreepBGRa;
    }

    double getTemperatureRelatedCoefficient(
        double const t, double const dt, ParameterLib::SpatialPosition const& x,
        double const T, double const deviatoric_stress_norm) const override;

private:
    NumLib::NewtonRaphsonSolverParameters const _nonlinear_solver_parameters;

    Parameter const& _a;        /// A parameter determined by experiment.
    Parameter const& _n;        /// Creep rate exponent n.
    Parameter const& _sigma_f;  /// A stress scaling factor.
    Parameter const& _q;        /// Activation energy
};

extern template class CreepBGRa<2>;
extern template class CreepBGRa<3>;

}  // end of namespace Creep
}  // end of namespace Solids
}  // namespace MaterialLib
