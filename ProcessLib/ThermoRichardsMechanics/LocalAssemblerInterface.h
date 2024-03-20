/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MaterialLib/SolidModels/SelectSolidConstitutiveRelation.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/Integration/GenericIntegrationMethod.h"
#include "ProcessLib/ConstitutiveRelations/StrainData.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "ProcessLib/Reflection/ReflectionSetIPData.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/Base.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/MaterialState.h"
#include "ProcessLib/ThermoRichardsMechanics/ThermoRichardsMechanicsProcessData.h"
#include "ProcessLib/Utils/SetOrGetIntegrationPointData.h"

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim, typename ConstitutiveTraits>
struct LocalAssemblerInterface : public ProcessLib::LocalAssemblerInterface,
                                 public NumLib::ExtrapolatableElement
{
    LocalAssemblerInterface(
        MeshLib::Element const& e,
        NumLib::GenericIntegrationMethod const& integration_method,
        bool const is_axially_symmetric,
        ThermoRichardsMechanicsProcessData<DisplacementDim, ConstitutiveTraits>&
            process_data)
        : process_data_(process_data),
          integration_method_(integration_method),
          element_(e),
          is_axially_symmetric_(is_axially_symmetric),
          solid_material_(MaterialLib::Solids::selectSolidConstitutiveRelation(
              process_data_.solid_materials, process_data_.material_ids,
              e.getID()))
    {
        unsigned const n_integration_points =
            integration_method_.getNumberOfPoints();

        current_states_.resize(n_integration_points);
        prev_states_.resize(n_integration_points);
        output_data_.resize(n_integration_points);

        material_states_.reserve(n_integration_points);
        for (unsigned ip = 0; ip < n_integration_points; ++ip)
        {
            material_states_.emplace_back(
                solid_material_.createMaterialStateVariables());

            // Set initial strain field to zero.
            std::get<StrainData<DisplacementDim>>(current_states_[ip]).eps =
                KelvinVector<DisplacementDim>::Zero();
        }
    }

    std::size_t setIPDataInitialConditions(std::string_view name,
                                           double const* values,
                                           int const integration_order)
    {
        if (integration_order !=
            static_cast<int>(integration_method_.getIntegrationOrder()))
        {
            OGS_FATAL(
                "Setting integration point initial conditions; The integration "
                "order of the local assembler for element {:d} is different "
                "from the integration order in the initial condition.",
                element_.getID());
        }

        if (name == "sigma" && process_data_.initial_stress.value)
        {
            OGS_FATAL(
                "Setting initial conditions for stress from integration "
                "point data and from a parameter '{:s}' is not possible "
                "simultaneously.",
                process_data_.initial_stress.value->name);
        }

        // TODO (naumov) this information is runtime information and I'm not
        // sure how to put it into the reflected data structure. The
        // reflectWithName function also supports only a single return value.
        if (name.starts_with("material_state_variable_"))
        {
            name.remove_prefix(24);

            auto const& internal_variables =
                solid_material_.getInternalVariables();
            if (auto const iv = std::find_if(
                    begin(internal_variables), end(internal_variables),
                    [&name](auto const& iv) { return iv.name == name; });
                iv != end(internal_variables))
            {
                DBUG("Setting material state variable '{:s}'", name);
                return ProcessLib::
                    setIntegrationPointDataMaterialStateVariables(
                        values, material_states_,
                        &MaterialStateData<
                            DisplacementDim>::material_state_variables,
                        iv->reference);
            }
            return 0;
        }

        // TODO this logic could be pulled out of the local assembler into the
        // process. That might lead to a slightly better performance due to less
        // string comparisons.
        return ProcessLib::Reflection::reflectSetIPData<DisplacementDim>(
            name, values, current_states_);
    }

    std::vector<double> getMaterialStateVariableInternalState(
        std::function<std::span<double>(
            typename MaterialLib::Solids::MechanicsBase<DisplacementDim>::
                MaterialStateVariables&)> const& get_values_span,
        int const& n_components) const
    {
        return ProcessLib::getIntegrationPointDataMaterialStateVariables(
            material_states_,
            &MaterialStateData<DisplacementDim>::material_state_variables,
            get_values_span, n_components);
    }

    // TODO move to NumLib::ExtrapolatableElement
    unsigned getNumberOfIntegrationPoints() const
    {
        return integration_method_.getNumberOfPoints();
    }

    int getMaterialID() const
    {
        return process_data_.material_ids == nullptr
                   ? 0
                   : (*process_data_.material_ids)[element_.getID()];
    }

    typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::MaterialStateVariables const&
    getMaterialStateVariablesAt(unsigned integration_point) const
    {
        return *material_states_[integration_point].material_state_variables;
    }

    void postTimestepConcrete(Eigen::VectorXd const& /*local_x*/,
                              Eigen::VectorXd const& /*local_x_prev*/,
                              double const /*t*/, double const /*dt*/,
                              int const /*process_id*/) override
    {
        unsigned const n_integration_points =
            integration_method_.getNumberOfPoints();

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            // TODO re-evaluate part of the assembly in order to be consistent?
            material_states_[ip].pushBackState();
        }

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            prev_states_[ip] = current_states_[ip];
        }
    }

    static auto getReflectionDataForOutput()
    {
        using Self =
            LocalAssemblerInterface<DisplacementDim, ConstitutiveTraits>;

        return ProcessLib::Reflection::reflectWithoutName(
            &Self::current_states_, &Self::output_data_);
    }

protected:
    ThermoRichardsMechanicsProcessData<DisplacementDim, ConstitutiveTraits>&
        process_data_;

    std::vector<typename ConstitutiveTraits::StatefulData>
        current_states_;  // TODO maybe do not store but rather re-evaluate for
                          // state update
    std::vector<typename ConstitutiveTraits::StatefulDataPrev> prev_states_;

    // Material state is special, because it contains both the current and the
    // old state.
    std::vector<MaterialStateData<DisplacementDim>> material_states_;

    NumLib::GenericIntegrationMethod const& integration_method_;
    MeshLib::Element const& element_;
    bool const is_axially_symmetric_;

    typename ConstitutiveTraits::SolidConstitutiveRelation const&
        solid_material_;

    std::vector<typename ConstitutiveTraits::OutputData> output_data_;
};

}  // namespace ProcessLib::ThermoRichardsMechanics
