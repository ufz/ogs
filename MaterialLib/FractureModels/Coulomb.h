/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <Eigen/Eigen>
#include <utility>

#include "NumLib/NewtonRaphson.h"
#include "ParameterLib/Parameter.h"

#include "FractureModelBase.h"

namespace MaterialLib
{
namespace Fracture
{
namespace Coulomb
{
template <int DisplacementDim>
struct StateVariables
    : public FractureModelBase<DisplacementDim>::MaterialStateVariables
{
    void setInitialConditions() { w_p = w_p_prev; }

    void pushBackState() override { w_p_prev = w_p; }

    /// Plastic component of the displacement jump in fracture's local
    /// coordinates.
    Eigen::Matrix<double, DisplacementDim, 1> w_p =
        Eigen::Matrix<double, DisplacementDim, 1>::Zero();

    // Initial values from previous timestep
    Eigen::Matrix<double, DisplacementDim, 1> w_p_prev;  ///< \copydoc w_p

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

template <int DisplacementDim>
class Coulomb final : public FractureModelBase<DisplacementDim>
{
public:
    std::unique_ptr<
        typename FractureModelBase<DisplacementDim>::MaterialStateVariables>
    createMaterialStateVariables() override
    {
        return std::make_unique<StateVariables<DisplacementDim>>();
    }

public:
    /// Variables specific to the material model
    struct MaterialProperties
    {
        using P = ParameterLib::Parameter<double>;
        using X = ParameterLib::SpatialPosition;

        MaterialProperties(
                P const& normal_stiffness_, P const& shear_stiffness_,
                P const& friction_angle_, P const& dilatancy_angle_,
                P const& cohesion_)
            : normal_stiffness(normal_stiffness_), shear_stiffness(shear_stiffness_),
              friction_angle(friction_angle_), dilatancy_angle(dilatancy_angle_),
              cohesion(cohesion_)
        {
        }

        /// Normal stiffness given in units of stress per length.
        P const& normal_stiffness;
        /// Shear stiffness given in units of stress per length.
        P const& shear_stiffness;
        /// Governs the normal stress dependence of the shear strength.
        /// \note Given in degrees (not radian).
        P const& friction_angle;
        /// Governs the dilatancy behaviour during the plastic deformation of
        /// the fault.
        /// \note Given in degrees (not radian).
        P const& dilatancy_angle;
        /// Fracture cohesion in units of stress.
        P const& cohesion;
    };


public:
    explicit Coulomb(
        NumLib::NewtonRaphsonSolverParameters nonlinear_solver_parameters,
        double const penalty_aperture_cutoff,
        bool const tension_cutoff,
        MaterialProperties material_properties)
        : _nonlinear_solver_parameters(std::move(nonlinear_solver_parameters)),
          _penalty_aperture_cutoff(penalty_aperture_cutoff),
          _tension_cutoff(tension_cutoff),
          _mp(std::move(material_properties))
    {
    }

    /**
     * Computation of the constitutive relation for the Mohr-Coulomb model.
     *
     * @param t           current time
     * @param x           current position in space
     * @param aperture0   initial fracture's aperture
     * @param sigma0      initial stress
     * @param w_prev      fracture displacement at previous time step
     * @param w           fracture displacement at current time step
     * @param sigma_prev  stress at previous time step
     * @param sigma       stress at current time step
     * @param Kep         tangent matrix for stress and fracture displacements
     * @param material_state_variables   material state variables
     */
    void computeConstitutiveRelation(
        double const t,
        ParameterLib::SpatialPosition const& x,
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
            Kep,
        typename FractureModelBase<DisplacementDim>::MaterialStateVariables&
            material_state_variables) override;

private:
    NumLib::NewtonRaphsonSolverParameters const _nonlinear_solver_parameters;

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

}  // namespace Coulomb
}  // namespace Fracture
}  // namespace MaterialLib

namespace MaterialLib
{
namespace Fracture
{
namespace Coulomb
{
extern template class Coulomb<2>;
extern template class Coulomb<3>;
}  // namespace Coulomb
}  // namespace Fracture
}  // namespace MaterialLib
