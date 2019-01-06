/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MaterialLib/SolidModels/MechanicsBase.h"

#include <MGIS/Behaviour/Behaviour.hxx>
#include <MGIS/Behaviour/BehaviourData.hxx>

#include "ProcessLib/Parameter/Parameter.h"

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
            : _data{b}
        {
        }

        MaterialStateVariables(MaterialStateVariables const&) = default;

        void pushBackState() override { mgis::behaviour::update(_data); }

        MaterialStateVariables& operator=(MaterialStateVariables const&) =
            default;

        typename MechanicsBase<DisplacementDim>::MaterialStateVariables&
        operator=(typename MechanicsBase<DisplacementDim>::
                      MaterialStateVariables const& state) noexcept override
        {
            return operator=(static_cast<MaterialStateVariables const&>(state));
        }

        // TODO fix: this has to be mutable.
        mutable mgis::behaviour::BehaviourData _data;
    };

    using KelvinVector =
        MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;
    using KelvinMatrix =
        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>;

    MFront(mgis::behaviour::Behaviour&& behaviour,
           std::vector<ProcessLib::Parameter<double> const*>&&
               material_properties);

    std::unique_ptr<
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables>
    createMaterialStateVariables() const override;

    boost::optional<std::tuple<KelvinVector,
                               std::unique_ptr<typename MechanicsBase<
                                   DisplacementDim>::MaterialStateVariables>,
                               KelvinMatrix>>
    integrateStress(
        double const t,
        ProcessLib::SpatialPosition const& x,
        double const dt,
        KelvinVector const& eps_prev,
        KelvinVector const& eps,
        KelvinVector const& sigma_prev,
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables const&
            material_state_variables,
        double const T) const override;

    double computeFreeEnergyDensity(
        double const t,
        ProcessLib::SpatialPosition const& x,
        double const dt,
        KelvinVector const& eps,
        KelvinVector const& sigma,
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables const&
            material_state_variables) const override;

private:
    mgis::behaviour::Behaviour _behaviour;
    std::vector<ProcessLib::Parameter<double> const*> _material_properties;
};

extern template class MFront<2>;
extern template class MFront<3>;

}  // namespace MFront
}  // namespace Solids
}  // namespace MaterialLib
