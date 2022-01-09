/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
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
            : _youngs_modulus(youngs_modulus), _poissons_ratio(poissons_ratio)
        {
        }

        /// Lamé's first parameter.
        double lambda(double const t, X const& x) const
        {
            return _youngs_modulus(t, x)[0] * _poissons_ratio(t, x)[0] /
                   (1 + _poissons_ratio(t, x)[0]) /
                   (1 - 2 * _poissons_ratio(t, x)[0]);
        }

        /// Lamé's second parameter, the shear modulus.
        double mu(double const t, X const& x) const
        {
            return _youngs_modulus(t, x)[0] /
                   (2 * (1 + _poissons_ratio(t, x)[0]));
        }

        /// the bulk modulus.
        double bulk_modulus(double const t, X const& x) const
        {
            return _youngs_modulus(t, x)[0] /
                   (3 * (1 - 2 * _poissons_ratio(t, x)[0]));
        }

    private:
        P const& _youngs_modulus;
        P const& _poissons_ratio;
    };

public:
    static int const KelvinVectorSize =
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);
    using KelvinVector =
        MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;
    using KelvinMatrix =
        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>;

    explicit LinearElasticIsotropic(MaterialProperties material_properties)
        : _mp(std::move(material_properties))
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
        MaterialPropertyLib::VariableArray const& variable_array_prev,
        MaterialPropertyLib::VariableArray const& variable_array,
        double const t, ParameterLib::SpatialPosition const& x,
        double const /*dt*/,
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables const&
            material_state_variables) const override;

    KelvinMatrix getElasticTensor(double const t,
                                  ParameterLib::SpatialPosition const& x,
                                  double const T) const;

    MaterialProperties getMaterialProperties() const { return _mp; }

    double getBulkModulus(double const t,
                          ParameterLib::SpatialPosition const& x,
                          KelvinMatrix const* const /*C*/) const override
    {
        return _mp.bulk_modulus(t, x);
    }



protected:
    MaterialProperties _mp;
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
