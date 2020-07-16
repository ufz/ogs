/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "MaterialLib/SolidModels/MechanicsBase.h"

#include <MGIS/Behaviour/Behaviour.hxx>
#include <MGIS/Behaviour/BehaviourData.hxx>

#include "ParameterLib/Parameter.h"

namespace MaterialLib
{
namespace Solids
{
namespace MFront
{
/// Converts MGIS kinematic to a string representation.
const char* toString(mgis::behaviour::Behaviour::Kinematic);

/// Converts MGIS symmetry to a string representation.
const char* toString(mgis::behaviour::Behaviour::Symmetry);

/// Converts MGIS behaviour type to a string representation.
const char* btypeToString(int);

/// Converts MGIS variable type to a string representation.
const char* varTypeToString(int);

/// Uses a material model provided by MFront (via MFront's generic interface and
/// the MGIS library).
template <int DisplacementDim>
class MFront : public MechanicsBase<DisplacementDim>
{
public:
    struct MaterialStateVariables
        : public MechanicsBase<DisplacementDim>::MaterialStateVariables
    {
        explicit MaterialStateVariables(mgis::behaviour::Behaviour const& b)
            : _behaviour_data{b}
        {
        }

        MaterialStateVariables(MaterialStateVariables const&) = default;
        MaterialStateVariables(MaterialStateVariables&&) = delete;

        void pushBackState() override
        {
            mgis::behaviour::update(_behaviour_data);
        }

        mgis::behaviour::BehaviourData _behaviour_data;
    };

    using KelvinVector =
        MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;
    using KelvinMatrix =
        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>;

    MFront(mgis::behaviour::Behaviour&& behaviour,
           std::vector<ParameterLib::Parameter<double> const*>&&
               material_properties,
           boost::optional<ParameterLib::CoordinateSystem> const&
               local_coordinate_system);

        std::unique_ptr<
            typename MechanicsBase<DisplacementDim>::
                MaterialStateVariables> createMaterialStateVariables()
            const override;

    std::optional<std::tuple<KelvinVector,
                             std::unique_ptr<typename MechanicsBase<
                                 DisplacementDim>::MaterialStateVariables>,
                             KelvinMatrix>>
    integrateStress(
        MaterialPropertyLib::VariableArray const& variable_array_prev,
        MaterialPropertyLib::VariableArray const& variable_array,
        double const t,
        ParameterLib::SpatialPosition const& x,
        double const dt,
        KelvinVector const& eps_prev,
        KelvinVector const& eps,
        KelvinVector const& sigma_prev,
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables const&
            material_state_variables,
        double const T) const override;

    std::vector<typename MechanicsBase<DisplacementDim>::InternalVariable>
    getInternalVariables() const override;

    double getBulkModulus(double const /*t*/,
                          ParameterLib::SpatialPosition const& /*x*/,
                          KelvinMatrix const* const /*C*/) const override;

    double computeFreeEnergyDensity(
        double const t,
        ParameterLib::SpatialPosition const& x,
        double const dt,
        KelvinVector const& eps,
        KelvinVector const& sigma,
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables const&
            material_state_variables) const override;

private:
    mgis::behaviour::Behaviour _behaviour;
    std::vector<ParameterLib::Parameter<double> const*> _material_properties;
    ParameterLib::CoordinateSystem const* const _local_coordinate_system;
};

extern template class MFront<2>;
extern template class MFront<3>;

}  // namespace MFront
}  // namespace Solids
}  // namespace MaterialLib
