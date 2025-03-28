/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <Eigen/Core>

#include "ParameterLib/Parameter.h"

namespace MaterialLib
{
namespace Fracture
{

/**
 * Interface for mechanical fracture material models. Provides updates of the
 * stress for a given current state and also a tangent at that position.
 */
template <int DisplacementDim>
class FractureModelBase
{
public:
    /// The MaterialStateVariables may store material model specific state
    /// (other than sigma and eps), which are usually material history
    /// dependent. The objects are stored by the user (usually in assembly per
    /// integration point) and are created via \ref
    /// createMaterialStateVariables().
    struct MaterialStateVariables
    {
        virtual ~MaterialStateVariables() = default;
        virtual void pushBackState() = 0;

        void reset()
        {
            _is_tensile_stress = false;
            _shear_yield_function_value = -1.0;
        }

        void setShearYieldFunctionValue(double Fs) { _shear_yield_function_value = Fs; }
        double getShearYieldFunctionValue() const { return _shear_yield_function_value; }

        void setTensileStress(bool flag) { _is_tensile_stress = flag; }
        bool setTensileStress() const { return _is_tensile_stress; }

    private:
        bool _is_tensile_stress = false;
        double _shear_yield_function_value = -1.0;
    };

    /// Polymorphic creator for MaterialStateVariables objects specific for a
    /// material model.
    virtual std::unique_ptr<MaterialStateVariables>
    createMaterialStateVariables() = 0;

    virtual ~FractureModelBase() = default;

    /**
     * Computation of the constitutive relation for specific material model.
     * This should be implemented in the derived model.
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
    virtual void computeConstitutiveRelation(
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
            C,
        MaterialStateVariables& material_state_variables) = 0;
};

}  // namespace Fracture
}  // namespace MaterialLib
