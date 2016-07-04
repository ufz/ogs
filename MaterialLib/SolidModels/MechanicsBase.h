/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MATERIALLIB_SOLIDMODELS_MECHANICSBASE_H_
#define MATERIALLIB_SOLIDMODELS_MECHANICSBASE_H_

#include <memory>

#include "ProcessLib/Deformation/BMatrixPolicy.h"
#include "ProcessLib/Parameter/Parameter.h"

namespace MeshLib
{
class Element;
}

namespace Solids
{
template <int DisplacementDim>
struct MechanicsBase
{
    /// The MaterialStateVariables may store material model specific state
    /// (other than sigma and eps), which are usually material history
    /// dependent. The objects are stored by the user (usually in assembly per
    /// integration point) and are created via \ref
    /// createMaterialStateVariables().
    struct MaterialStateVariables
    {
        virtual ~MaterialStateVariables() = default;
        virtual void pushBackState() = 0;
    };

    using KelvinVector = ProcessLib::KelvinVectorType<DisplacementDim>;
    using KelvinMatrix = ProcessLib::KelvinMatrixType<DisplacementDim>;

    void computeConstitutiveRelation(
        double const t,
        ProcessLib::SpatialPosition const& x,
        double const dt,
        Eigen::Matrix<double, Eigen::Dynamic, 1> const& eps_prev,
        Eigen::Matrix<double, Eigen::Dynamic, 1> const& eps,
        Eigen::Matrix<double, Eigen::Dynamic, 1> const& sigma_prev,
        Eigen::Matrix<double, Eigen::Dynamic, 1>& sigma,
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>&
            C,
        MaterialStateVariables& material_state_variables)
    {
        // TODO Avoid copies of data:
        // Using MatrixBase<Derived> not possible because template functions
        // cannot be virtual. Maybe there is a workaround for this.  Using
        // Map<Matrix<double, ...>> makes the interface (for the material model
        // implementation) unnecessary difficult.
        KelvinVector const eps_prev_{eps_prev};
        KelvinVector const eps_{eps};
        KelvinVector const sigma_prev_{sigma_prev};
        KelvinVector sigma_{sigma};
        KelvinMatrix C_{C};

        computeConstitutiveRelation(t,
                                    x,
                                    dt,
                                    eps_prev_,
                                    eps_,
                                    sigma_prev_,
                                    sigma_,
                                    C_,
                                    material_state_variables);

        sigma = sigma_;
        C = C_;
    }

    virtual void computeConstitutiveRelation(
        double const t,
        ProcessLib::SpatialPosition const& x,
        double const dt,
        KelvinVector const& eps_prev,
        KelvinVector const& eps,
        KelvinVector const& sigma_prev,
        KelvinVector& sigma,
        KelvinMatrix& C,
        MaterialStateVariables& material_state_variables) = 0;

    virtual std::unique_ptr<MaterialStateVariables>
    createMaterialStateVariables() = 0;

    virtual ~MechanicsBase() = default;
};
}

#endif  // MATERIALLIB_SOLIDMODELS_MECHANICSBASE_H_
