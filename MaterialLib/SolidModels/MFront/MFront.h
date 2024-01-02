/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "MFrontGeneric.h"
#include "Variable.h"

namespace MaterialLib
{
namespace Solids
{
namespace MFront
{
template <int DisplacementDim>
class MFront : private MFrontGeneric<DisplacementDim,
                                     boost::mp11::mp_list<Strain>,
                                     boost::mp11::mp_list<Stress>,
                                     boost::mp11::mp_list<Temperature>>,
               public MechanicsBase<DisplacementDim>
{
    using Base = MFrontGeneric<DisplacementDim,
                               boost::mp11::mp_list<Strain>,
                               boost::mp11::mp_list<Stress>,
                               boost::mp11::mp_list<Temperature>>;
    using KelvinVector = typename Base::KelvinVector;
    using KelvinMatrix = typename Base::KelvinMatrix;

public:
    using Base::Base;

    std::unique_ptr<
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables>
    createMaterialStateVariables() const override
    {
        return Base::createMaterialStateVariables();
    }

    void initializeInternalStateVariables(
        double const t,
        ParameterLib::SpatialPosition const& x,
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables&
            material_state_variables) const override
    {
        Base::initializeInternalStateVariables(t, x, material_state_variables);
    }

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
            material_state_variables) const override
    {
        auto res = Base::integrateStress(variable_array_prev,
                                         variable_array,
                                         t,
                                         x,
                                         dt,
                                         material_state_variables);

        if (!res)
        {
            return std::nullopt;
        }

        auto& [stress_data, state, tangent_operator_data] = *res;
        auto const view = this->createThermodynamicForcesView();

        auto const C =
            blocks_view_.block(stress, strain, tangent_operator_data);

        return std::optional<
            std::tuple<KelvinVector,
                       std::unique_ptr<typename MechanicsBase<
                           DisplacementDim>::MaterialStateVariables>,
                       KelvinMatrix>>{std::in_place,
                                      view.block(stress, stress_data),
                                      std::move(state), C};
    }

    std::vector<typename MechanicsBase<DisplacementDim>::InternalVariable>
    getInternalVariables() const override
    {
        return Base ::getInternalVariables();
    }

    double getBulkModulus(double const t,
                          ParameterLib::SpatialPosition const& x,
                          KelvinMatrix const* const C) const override
    {
        return Base::getBulkModulus(t, x, C);
    }

    double computeFreeEnergyDensity(
        double const t,
        ParameterLib::SpatialPosition const& x,
        double const dt,
        KelvinVector const& eps,
        KelvinVector const& sigma,
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables const&
            material_state_variables) const override
    {
        return Base::computeFreeEnergyDensity(
            t, x, dt, eps, sigma, material_state_variables);
    }

private:
    OGSMFrontTangentOperatorBlocksView<
        DisplacementDim,
        ForcesGradsCombinations<boost::mp11::mp_list<Strain>,
                                boost::mp11::mp_list<Stress>,
                                boost::mp11::mp_list<Temperature>>::type>
        blocks_view_ = this->createTangentOperatorBlocksView();
};

extern template class MFront<2>;
extern template class MFront<3>;

}  // namespace MFront
}  // namespace Solids
}  // namespace MaterialLib
