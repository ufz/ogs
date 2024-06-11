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

#include "ConstitutiveRelations/ConstitutiveData.h"
#include "MaterialLib/SolidModels/MechanicsBase.h"
#include "MaterialLib/SolidModels/SelectSolidConstitutiveRelation.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/Integration/GenericIntegrationMethod.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "ProcessLib/Reflection/ReflectionSetIPData.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/MaterialState.h"
#include "ProcessLib/Utils/SetOrGetIntegrationPointData.h"
#include "RichardsMechanicsProcessData.h"

namespace ProcessLib
{
namespace RichardsMechanics
{
template <int DisplacementDim>
struct LocalAssemblerInterface : public ProcessLib::LocalAssemblerInterface,
                                 public NumLib::ExtrapolatableElement
{
    LocalAssemblerInterface(
        MeshLib::Element const& e,
        NumLib::GenericIntegrationMethod const& integration_method,
        bool const is_axially_symmetric,
        RichardsMechanicsProcessData<DisplacementDim>& process_data)
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

        if (name == "sigma" && process_data_.initial_stress != nullptr)
        {
            OGS_FATAL(
                "Setting initial conditions for stress from integration "
                "point data and from a parameter '{:s}' is not possible "
                "simultaneously.",
                process_data_.initial_stress->name);
        }

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
                        &ProcessLib::ThermoRichardsMechanics::MaterialStateData<
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

    virtual std::vector<double> getMaterialStateVariableInternalState(
        std::function<std::span<double>(
            typename MaterialLib::Solids::MechanicsBase<DisplacementDim>::
                MaterialStateVariables&)> const& get_values_span,
        int const& n_components) const = 0;

    // TODO move to NumLib::ExtrapolatableElement
    virtual unsigned getNumberOfIntegrationPoints() const = 0;

    virtual int getMaterialID() const = 0;

    virtual typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::MaterialStateVariables const&
    getMaterialStateVariablesAt(unsigned /*integration_point*/) const = 0;

    static auto getReflectionDataForOutput()
    {
        using Self = LocalAssemblerInterface<DisplacementDim>;

        return ProcessLib::Reflection::reflectWithoutName(
            &Self::current_states_, &Self::output_data_);
    }

protected:
    RichardsMechanicsProcessData<DisplacementDim>& process_data_;
    NumLib::GenericIntegrationMethod const& integration_method_;
    MeshLib::Element const& element_;
    bool const is_axially_symmetric_;

    MaterialLib::Solids::MechanicsBase<DisplacementDim> const& solid_material_;

    std::vector<StatefulData<DisplacementDim>> current_states_;
    std::vector<StatefulDataPrev<DisplacementDim>> prev_states_;
    std::vector<
        ProcessLib::ThermoRichardsMechanics::MaterialStateData<DisplacementDim>>
        material_states_;
    std::vector<OutputData<DisplacementDim>> output_data_;
};
}  // namespace RichardsMechanics
}  // namespace ProcessLib
