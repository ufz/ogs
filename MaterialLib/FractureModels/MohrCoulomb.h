/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/Eigen>

#include "ProcessLib/Parameter/Parameter.h"

#include "FractureModelBase.h"

namespace MaterialLib
{
namespace Fracture
{

/**
 * Mohr-Coulomb fracture model
 */
template <int DisplacementDim>
class MohrCoulomb final : public FractureModelBase<DisplacementDim>
{
public:
    /// Variables specific to the material model
    struct MaterialProperties
    {
        using P = ProcessLib::Parameter<double>;
        using X = ProcessLib::SpatialPosition;

        MaterialProperties(
                P const& normal_stiffness, P const& shear_stiffness,
                P const& friction_angle, P const& dilatancy_angle,
                P const& cohesion)
            : normal_stiffness(normal_stiffness), shear_stiffness(shear_stiffness),
              friction_angle(friction_angle), dilatancy_angle(dilatancy_angle),
              cohesion(cohesion)
        {
        }

        P const& normal_stiffness;
        P const& shear_stiffness;
        P const& friction_angle;
        P const& dilatancy_angle;
        P const& cohesion;
    };

public:

    explicit MohrCoulomb(
        MaterialProperties const& material_properties)
        : _mp(material_properties)
    {
    }

    void computeConstitutiveRelation(
            double const t,
            ProcessLib::SpatialPosition const& x,
            Eigen::Ref<Eigen::VectorXd const> w_prev,
            Eigen::Ref<Eigen::VectorXd const> w,
            Eigen::Ref<Eigen::VectorXd const> sigma_prev,
            Eigen::Ref<Eigen::VectorXd> sigma,
            Eigen::Ref<Eigen::MatrixXd> Kep) override;

private:

    MaterialProperties _mp;
};

}  // namespace Fracture
}  // namespace MaterialLib

namespace MaterialLib
{
namespace Fracture
{
extern template class MohrCoulomb<2>;
extern template class MohrCoulomb<3>;
}  // namespace Fractrue
}  // namespace MaterialLib

