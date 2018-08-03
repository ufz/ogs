/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file   CreepBGRa.h
 *  Created on July 6, 2018, 9:53 AM
 */

#pragma once

#include <cmath>
#include <memory>
#include <tuple>

#include "LinearElasticIsotropic.h"
#include "MathLib/KelvinVector.h"
#include "NumLib/NewtonRaphson.h"

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
 *
 */
template <int DisplacementDim>
class CreepBGRa final : public LinearElasticIsotropic<DisplacementDim>
{
public:
    static int const KelvinVectorSize =
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value;
    using ResidualVectorType = Eigen::Matrix<double, KelvinVectorSize, 1>;
    using JacobianMatrix = Eigen::Matrix<double, KelvinVectorSize,
                                         KelvinVectorSize, Eigen::RowMajor>;

    using KelvinVector =
        MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;
    using KelvinMatrix =
        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>;

    std::unique_ptr<
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables>
    createMaterialStateVariables() override
    {
        return LinearElasticIsotropic<
            DisplacementDim>::createMaterialStateVariables();
    }

    CreepBGRa(
        typename LinearElasticIsotropic<DisplacementDim>::MaterialProperties mp,
        NumLib::NewtonRaphsonSolverParameters nonlinear_solver_parameters,
        const double A, const double n, const double sigma0, const double Q)
        : LinearElasticIsotropic<DisplacementDim>(std::move(mp)),
          _nonlinear_solver_parameters(std::move(nonlinear_solver_parameters)),
          _n(n),
          _coef(A * std::pow(1.5, 0.5 * (1 + _n)) / std::pow(sigma0, _n)),
          _q(Q)
    {
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
        double const T) override;

    ConstitutiveModel getConstitutiveModel() const override
    {
        return ConstitutiveModel::CreepBGRa;
    }

    double getTemperatureRelatedCoefficient(
        double const t, double const dt, ProcessLib::SpatialPosition const& x,
        double const T, double const deviatoric_stress_norm) const override;

private:
    NumLib::NewtonRaphsonSolverParameters const _nonlinear_solver_parameters;

    const double _n;  /// Creep rate exponent n.
    /// \f$A\left(\frac{3}{2}\right)^{n/2+1}/\sigma_{eff}^n \f$
    const double _coef;
    const double _q;  /// Activation energy
};

extern template class CreepBGRa<2>;
extern template class CreepBGRa<3>;

}  // end of namespace Creep
}  // end of namespace Solids
}  // namespace MaterialLib
