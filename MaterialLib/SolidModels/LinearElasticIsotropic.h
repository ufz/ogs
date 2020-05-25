/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "MechanicsBase.h"
#include "ParameterLib/Parameter.h"

namespace MaterialLib
{
namespace Solids
{
template <int DisplacementDim>
class LinearElasticIsotropic : public MechanicsBase<DisplacementDim>
{
public:
    /// Variables specific to the material model
    class MaterialProperties
    {
        using P = ParameterLib::Parameter<double>;
        using X = ParameterLib::SpatialPosition;

    public:
        MaterialProperties(P const& youngs_modulus, P const& poissons_ratio)
            : youngs_modulus_(youngs_modulus), poissons_ratio_(poissons_ratio)
        {
        }

        /// Lamé's first parameter.
        double lambda(double const t, X const& x) const
        {
            return youngs_modulus_(t, x)[0] * poissons_ratio_(t, x)[0] /
                   (1 + poissons_ratio_(t, x)[0]) /
                   (1 - 2 * poissons_ratio_(t, x)[0]);
        }

        /// Lamé's second parameter, the shear modulus.
        double mu(double const t, X const& x) const
        {
            return youngs_modulus_(t, x)[0] /
                   (2 * (1 + poissons_ratio_(t, x)[0]));
        }

        /// the bulk modulus.
        double bulk_modulus(double const t, X const& x) const
        {
            return youngs_modulus_(t, x)[0] /
                   (3 * (1 - 2 * poissons_ratio_(t, x)[0]));
        }

    private:
        P const& youngs_modulus_;
        P const& poissons_ratio_;
    };

public:
    static int const KelvinVectorSize =
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value;
    using KelvinVector =
        MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;
    using KelvinMatrix =
        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>;

    explicit LinearElasticIsotropic(MaterialProperties material_properties)
        : mp_(std::move(material_properties))
    {
    }

    double computeFreeEnergyDensity(
        double const /*t*/,
        ParameterLib::SpatialPosition const& /*x*/,
        double const /*dt*/,
        KelvinVector const& eps,
        KelvinVector const& sigma,
        typename MechanicsBase<DisplacementDim>::
            MaterialStateVariables const& /* material_state_variables */)
        const override
    {
        return eps.dot(sigma) / 2;
    }

    std::optional<
        std::tuple<typename MechanicsBase<DisplacementDim>::KelvinVector,
                   std::unique_ptr<typename MechanicsBase<
                       DisplacementDim>::MaterialStateVariables>,
                   typename MechanicsBase<DisplacementDim>::KelvinMatrix>>
    integrateStress(
        double const t, ParameterLib::SpatialPosition const& x,
        double const /*dt*/, KelvinVector const& eps_prev,
        KelvinVector const& eps, KelvinVector const& sigma_prev,
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables const&
            material_state_variables,
        double const T) const override;

    KelvinMatrix getElasticTensor(double const t,
                                  ParameterLib::SpatialPosition const& x,
                                  double const T) const;

    MaterialProperties getMaterialProperties() const { return mp_; }

    double getBulkModulus(double const t,
                          ParameterLib::SpatialPosition const& x,
                          KelvinMatrix const* const /*C*/) const override
    {
        return mp_.bulk_modulus(t, x);
    }



protected:
    MaterialProperties mp_;
};

extern template class LinearElasticIsotropic<2>;
extern template class LinearElasticIsotropic<3>;

template <int DisplacementDim>
MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>
elasticTangentStiffness(double const first_lame_parameter,
                        double const shear_modulus)
{
    using KelvinMatrix =
        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>;

    KelvinMatrix tangentStiffness = KelvinMatrix::Zero();
    tangentStiffness.template topLeftCorner<3, 3>().setConstant(
        first_lame_parameter);
    tangentStiffness.noalias() += 2 * shear_modulus * KelvinMatrix::Identity();
    return tangentStiffness;
}

}  // namespace Solids
}  // namespace MaterialLib
