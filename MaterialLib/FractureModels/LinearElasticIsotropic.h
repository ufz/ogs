/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
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

        /// Normal stiffness given in units of stress per length.
        P const& normal_stiffness;
        /// Shear stiffness given in units of stress per length.
        P const& shear_stiffness;
    };

    struct MaterialStateVariables
        : public FractureModelBase<DisplacementDim>::MaterialStateVariables
    {
        void pushBackState() override {}
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
    explicit LinearElasticIsotropic(double const penalty_aperture_cutoff,
                                    bool const tension_cutoff,
                                    MaterialProperties material_properties)
        : _penalty_aperture_cutoff(penalty_aperture_cutoff),
          _tension_cutoff(tension_cutoff),
          _mp(std::move(material_properties))
    {
    }

    /**
     * Computation of the constitutive relation for the linear elastic model.
     *
     * @param t           current time
     * @param x           current position in space
     * @param aperture0   initial fracture's aperture
     * @param sigma0      initial stress
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
        double const aperture0,
        Eigen::Ref<Eigen::VectorXd const>
            sigma0,
        Eigen::Ref<Eigen::VectorXd const>
            w_prev,
        Eigen::Ref<Eigen::VectorXd const>
            w,
        Eigen::Ref<Eigen::VectorXd const>
            sigma_prev,
        Eigen::Ref<Eigen::VectorXd>
            sigma,
        Eigen::Ref<Eigen::MatrixXd>
            C,
        typename FractureModelBase<DisplacementDim>::MaterialStateVariables&
            material_state_variables) override;

private:
    /// Compressive normal displacements above this value will not enter the
    /// computation of the normal stiffness modulus of the fracture.
    /// \note Setting this to the initial aperture value allows negative
    /// apertures.
    double const _penalty_aperture_cutoff;

    /// If set no resistance to open the fracture over the initial aperture is
    /// opposed.
    bool const _tension_cutoff;

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
}  // namespace Fracture
}  // namespace MaterialLib
