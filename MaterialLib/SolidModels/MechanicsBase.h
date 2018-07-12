/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <boost/optional.hpp>
#include <functional>
#include <memory>
#include <tuple>
#include <vector>

#include "BaseLib/Error.h"
#include "ProcessLib/Deformation/BMatrixPolicy.h"

namespace ProcessLib
{
class SpatialPosition;
}

namespace MeshLib
{
class Element;
}

namespace MaterialLib
{
namespace Solids
{
enum class ConstitutiveModel
{
    Ehlers,
    LinearElasticIsotropic,
    Lubby2,
    CreepBGRa,
    Invalid
};

/// Interface for mechanical solid material models. Provides updates of the
/// stress for a given current state and also a tangent at that position. If the
/// implemented material model stores an internal state, the nested
/// MaterialStateVariables class should be used; it's only responsibility is to
/// provide state's push back possibility.
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
        virtual MaterialStateVariables& operator=(
            MaterialStateVariables const&) = default;

        virtual void pushBackState() = 0;
    };

    /// Polymorphic creator for MaterialStateVariables objects specific for a
    /// material model.
    virtual std::unique_ptr<MaterialStateVariables>
    createMaterialStateVariables() = 0;

    using KelvinVector =
        MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;
    using KelvinMatrix =
        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>;

    /// Dynamic size Kelvin vector and matrix wrapper for the polymorphic
    /// constitutive relation compute function.
    /// Returns nothing in case of errors in the computation if Newton
    /// iterations did not converge, for example.
    boost::optional<std::tuple<KelvinVector,
                               std::unique_ptr<MaterialStateVariables>,
                               KelvinMatrix>>
    integrateStress(double const t,
                    ProcessLib::SpatialPosition const& x,
                    double const dt,
                    Eigen::Matrix<double, Eigen::Dynamic, 1> const& eps_prev,
                    Eigen::Matrix<double, Eigen::Dynamic, 1> const& eps,
                    Eigen::Matrix<double, Eigen::Dynamic, 1> const& sigma_prev,
                    MaterialStateVariables const& material_state_variables,
                    double const T)
    {
        // TODO Avoid copies of data:
        // Using MatrixBase<Derived> not possible because template functions
        // cannot be virtual. Maybe there is a workaround for this.  Using
        // Map<Matrix<double, ...>> makes the interface (for the material model
        // implementation) unnecessary difficult.
        KelvinVector const eps_prev_{eps_prev};
        KelvinVector const eps_{eps};
        KelvinVector const sigma_prev_{sigma_prev};

        return integrateStress(
            t, x, dt, eps_prev_, eps_, sigma_prev_, material_state_variables,
            T);
    }

    /// Computation of the constitutive relation for specific material model.
    /// This should be implemented in the derived model. Fixed Kelvin vector and
    /// matrix size version; for dynamic size arguments there is an overloaded
    /// wrapper function.
    /// Returns nothing in case of errors in the computation if Newton
    /// iterations did not converge, for example.
    virtual boost::optional<std::tuple<KelvinVector,
                                       std::unique_ptr<MaterialStateVariables>,
                                       KelvinMatrix>>
    integrateStress(double const t,
                    ProcessLib::SpatialPosition const& x,
                    double const dt,
                    KelvinVector const& eps_prev,
                    KelvinVector const& eps,
                    KelvinVector const& sigma_prev,
                    MaterialStateVariables const& material_state_variables,
                    double const T) = 0;

    /// Helper type for providing access to internal variables.
    struct InternalVariable
    {
        using Getter = std::function<std::vector<double> const&(
            MaterialStateVariables const&, std::vector<double>& /*cache*/)>;

        /// name of the internal variable
        std::string const name;

        /// number of components of the internal variable
        unsigned const num_components;

        /// function accessing the internal variable
        Getter const getter;
    };

    /// Returns internal variables defined by the specific material model, if
    /// any.
    virtual std::vector<InternalVariable> getInternalVariables() const
    {
        return {};
    }

    /// Gets the type of constitutive model
    virtual ConstitutiveModel getConstitutiveModel() const
    {
        return ConstitutiveModel::Invalid;
    }

    /// Get temperature related coefficient for the global assembly if there is
    /// one.
    virtual double getTemperatureRelatedCoefficient(
        double const /*t*/, double const /*dt*/,
        ProcessLib::SpatialPosition const& /*x*/, double const /*T*/,
        double const /*deviatoric_stress_norm*/) const
    {
        return 0.0;
    }

    virtual double computeFreeEnergyDensity(
        double const t,
        ProcessLib::SpatialPosition const& x,
        double const dt,
        KelvinVector const& eps,
        KelvinVector const& sigma,
        MaterialStateVariables const& material_state_variables) const = 0;

    virtual ~MechanicsBase() = default;
};

}  // namespace Solids
}  // namespace MaterialLib
