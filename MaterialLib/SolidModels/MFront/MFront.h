/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
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
        explicit MaterialStateVariables(
            int const equivalent_plastic_strain_offset,
            mgis::behaviour::Behaviour const& b)
            : equivalent_plastic_strain_offset_(
                  equivalent_plastic_strain_offset),
              _behaviour_data{b}
        {
        }

        MaterialStateVariables(MaterialStateVariables const&) = default;
        MaterialStateVariables(MaterialStateVariables&&) = delete;

        void pushBackState() override
        {
            mgis::behaviour::update(_behaviour_data);
        }

        int const equivalent_plastic_strain_offset_;
        mgis::behaviour::BehaviourData _behaviour_data;

        double getEquivalentPlasticStrain() const override;
    };

    using KelvinVector =
        MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;
    using KelvinMatrix =
        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>;

    MFront(mgis::behaviour::Behaviour&& behaviour,
           std::vector<ParameterLib::Parameter<double> const*>&&
               material_properties,
           std::optional<ParameterLib::CoordinateSystem> const&
               local_coordinate_system);

    std::unique_ptr<
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables>
    createMaterialStateVariables() const override;

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
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables const&
            material_state_variables) const override;

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
    int const equivalent_plastic_strain_offset_;
    std::vector<ParameterLib::Parameter<double> const*> _material_properties;
    ParameterLib::CoordinateSystem const* const _local_coordinate_system;
};

extern template class MFront<2>;
extern template class MFront<3>;

}  // namespace MFront
}  // namespace Solids
}  // namespace MaterialLib
