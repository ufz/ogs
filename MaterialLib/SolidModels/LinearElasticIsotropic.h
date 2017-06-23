/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MechanicsBase.h"
#include "ProcessLib/Parameter/Parameter.h"

namespace MaterialLib
{
namespace Solids
{
template <int DisplacementDim>
class LinearElasticIsotropic final : public MechanicsBase<DisplacementDim>
{
public:
    /// Variables specific to the material model
    class MaterialProperties
    {
        using P = ProcessLib::Parameter<double>;
        using X = ProcessLib::SpatialPosition;

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

    private:
        P const& _youngs_modulus;
        P const& _poissons_ratio;
    };

    struct MaterialStateVariables
        : public MechanicsBase<DisplacementDim>::MaterialStateVariables
    {
        void pushBackState() override {}

        MaterialStateVariables& operator=(MaterialStateVariables const&) =
            default;
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables&
        operator=(typename MechanicsBase<DisplacementDim>::
                      MaterialStateVariables const& state) noexcept override
        {
            return operator=(static_cast<MaterialStateVariables const&>(state));
        }
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
    static int const KelvinVectorSize =
        ProcessLib::KelvinVectorDimensions<DisplacementDim>::value;
    using KelvinVector = ProcessLib::KelvinVectorType<DisplacementDim>;
    using KelvinMatrix = ProcessLib::KelvinMatrixType<DisplacementDim>;

    explicit LinearElasticIsotropic(MaterialProperties material_properties)
        : _mp(std::move(material_properties))
    {
    }

    boost::optional<std::tuple<KelvinVector,
                               std::unique_ptr<typename MechanicsBase<
                                   DisplacementDim>::MaterialStateVariables>,
                               KelvinMatrix>>
    integrateStress(
        double const t,
        ProcessLib::SpatialPosition const& x,
        double const /*dt*/,
        KelvinVector const& eps_prev,
        KelvinVector const& eps,
        KelvinVector const& sigma_prev,
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables const&
            material_state_variables) override;

private:
    MaterialProperties _mp;
};

extern template class LinearElasticIsotropic<2>;
extern template class LinearElasticIsotropic<3>;

}  // namespace Solids
}  // namespace MaterialLib
