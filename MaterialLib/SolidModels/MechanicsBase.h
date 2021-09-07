/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <functional>
#include <memory>
#include <optional>
#include <tuple>
#include <vector>

#include "BaseLib/DynamicSpan.h"
#include "BaseLib/Error.h"
#include "MaterialLib/MPL/VariableType.h"
#include "MathLib/KelvinVector.h"

namespace ParameterLib
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
        virtual void pushBackState(){};
        virtual double getEquivalentPlasticStrain() const { return 0.0; }
    };

    /// Polymorphic creator for MaterialStateVariables objects specific for a
    /// material model.
    virtual std::unique_ptr<MaterialStateVariables>
    createMaterialStateVariables() const
    {
        return std::make_unique<MaterialStateVariables>();
    }

    using KelvinVector =
        MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;
    using KelvinMatrix =
        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>;

    /// Computation of the constitutive relation for specific material model.
    /// This should be implemented in the derived model. Fixed Kelvin vector and
    /// matrix size version; for dynamic size arguments there is an overloaded
    /// wrapper function.
    /// Returns nothing in case of errors in the computation if Newton
    /// iterations did not converge, for example.
    virtual std::optional<std::tuple<
        KelvinVector, std::unique_ptr<MaterialStateVariables>, KelvinMatrix>>
    integrateStress(
        MaterialPropertyLib::VariableArray const& variable_array_prev,
        MaterialPropertyLib::VariableArray const& variable_array,
        double const t,
        ParameterLib::SpatialPosition const& x,
        double const dt,
        MaterialStateVariables const& material_state_variables) const = 0;

    /// Helper type for providing access to internal variables.
    struct InternalVariable
    {
        using Getter = std::function<std::vector<double> const&(
            MaterialStateVariables const&, std::vector<double>& /*cache*/)>;
        using WriteAccess = std::function<BaseLib::DynamicSpan<double>(
            MaterialStateVariables&)>;

        /// name of the internal variable
        std::string const name;

        /// number of components of the internal variable
        int const num_components;

        /// function accessing the internal variable
        Getter const getter;

        /// function accessing the internal variable
        WriteAccess const reference;
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

    virtual double getBulkModulus(
        double const /*t*/,
        ParameterLib::SpatialPosition const& /*x*/,
        KelvinMatrix const* const /*C*/ = nullptr) const
    {
        OGS_FATAL(
            "getBulkModulus is not yet implemented for this Solid Material.");
    }

    /// Get temperature related coefficient for the global assembly if there is
    /// one.
    virtual double getTemperatureRelatedCoefficient(
        double const /*t*/, double const /*dt*/,
        ParameterLib::SpatialPosition const& /*x*/, double const /*T*/,
        double const /*deviatoric_stress_norm*/) const
    {
        return 0.0;
    }

    virtual double computeFreeEnergyDensity(
        double const t,
        ParameterLib::SpatialPosition const& x,
        double const dt,
        KelvinVector const& eps,
        KelvinVector const& sigma,
        MaterialStateVariables const& material_state_variables) const = 0;

    virtual ~MechanicsBase() = default;
};

}  // namespace Solids
}  // namespace MaterialLib
