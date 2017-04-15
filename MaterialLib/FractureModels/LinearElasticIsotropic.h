/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/Eigen>
#include <utility>

#include "ProcessLib/Parameter/Parameter.h"

#include "FractureModelBase.h"

namespace MaterialLib
{
namespace Fracture
{

template <int DisplacementDim>
class LinearElasticIsotropic final : public FractureModelBase<DisplacementDim>
{
public:
    /// Variables specific to the material model
    struct MaterialProperties
    {
        using P = ProcessLib::Parameter<double>;
        using X = ProcessLib::SpatialPosition;

        MaterialProperties(P const& normal_stiffness_, P const& shear_stiffness_)
            : normal_stiffness(normal_stiffness_), shear_stiffness(shear_stiffness_)
        {
        }

        P const& normal_stiffness;
        P const& shear_stiffness;
    };

    struct MaterialStateVariables
        : public FractureModelBase<DisplacementDim>::MaterialStateVariables
    {
        void pushBackState() {}
    };

    std::unique_ptr<
        typename FractureModelBase<DisplacementDim>::MaterialStateVariables>
    createMaterialStateVariables() override
    {
        return std::unique_ptr<typename FractureModelBase<
            DisplacementDim>::MaterialStateVariables>{
            new MaterialStateVariables};
    }

public:
    explicit LinearElasticIsotropic(MaterialProperties material_properties)
        : _mp(std::move(material_properties))
    {
    }

    /**
     * Computation of the constitutive relation for the linear elastic model.
     *
     * @param t           current time
     * @param x           current position in space
     * @param w_prev      fracture displacement at previous time step
     * @param w           fracture displacement at current time step
     * @param sigma_prev  stress at previous time step
     * @param sigma       stress at current time step
     * @param C           tangent matrix for stress and fracture displacements
     * @param material_state_variables   material state variables
     */
    void computeConstitutiveRelation(
            double const t,
            ProcessLib::SpatialPosition const& x,
            Eigen::Ref<Eigen::VectorXd const> w_prev,
            Eigen::Ref<Eigen::VectorXd const> w,
            Eigen::Ref<Eigen::VectorXd const> sigma_prev,
            Eigen::Ref<Eigen::VectorXd> sigma,
            Eigen::Ref<Eigen::MatrixXd> C,
            typename FractureModelBase<DisplacementDim>::MaterialStateVariables&
            material_state_variables) override;

private:
    MaterialProperties _mp;
};

}  // namespace Fracture
}  // namespace MaterialLib

namespace MaterialLib
{
namespace Fracture
{
extern template class LinearElasticIsotropic<2>;
extern template class LinearElasticIsotropic<3>;
}  // namespace Fractrue
}  // namespace MaterialLib
