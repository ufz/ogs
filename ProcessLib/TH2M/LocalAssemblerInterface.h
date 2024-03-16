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
#include "ConstitutiveRelations/MaterialState.h"
#include "ConstitutiveRelations/SolidMechanics.h"
#include "MaterialLib/SolidModels/SelectSolidConstitutiveRelation.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/Integration/GenericIntegrationMethod.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "ProcessLib/Reflection/ReflectionSetIPData.h"
#include "ProcessLib/Utils/SetOrGetIntegrationPointData.h"
#include "TH2MProcessData.h"

namespace ProcessLib
{
namespace TH2M
{
template <int DisplacementDim>
struct LocalAssemblerInterface : public ProcessLib::LocalAssemblerInterface,
                                 public NumLib::ExtrapolatableElement
{
    LocalAssemblerInterface(
        MeshLib::Element const& e,
        NumLib::GenericIntegrationMethod const& integration_method,
        bool const is_axially_symmetric,
        TH2MProcessData<DisplacementDim>& process_data)
        : process_data_(process_data),
          integration_method_(integration_method),
          element_(e),
          is_axially_symmetric_(is_axially_symmetric),
          solid_material_(MaterialLib::Solids::selectSolidConstitutiveRelation(
              process_data_.solid_materials, process_data_.material_ids,
              element_.getID()))
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

            // Set initial stress field to zero. Might be overwritten by
            // integration point data or initial stress.
            current_states_[ip].eff_stress_data.sigma =
                KelvinVector<DisplacementDim>::Zero();
        }
    }

    virtual std::size_t setIPDataInitialConditions(
        std::string_view name, double const* values,
        int const integration_order) = 0;

    virtual std::vector<double> const& getIntPtDarcyVelocityGas(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtDarcyVelocityLiquid(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtDiffusionVelocityVapourGas(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtDiffusionVelocityGasGas(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtDiffusionVelocitySoluteLiquid(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtDiffusionVelocityLiquidLiquid(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtLiquidDensity(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtGasDensity(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtSolidDensity(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtVapourPressure(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtPorosity(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtMoleFractionGas(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtMassFractionGas(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtMassFractionLiquid(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtEnthalpyGas(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtEnthalpyLiquid(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtEnthalpySolid(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

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

    std::vector<double> getMaterialStateVariableInternalState(
        std::function<std::span<double>(
            typename ConstitutiveRelations::SolidConstitutiveRelation<
                DisplacementDim>::MaterialStateVariables&)> const&
            get_values_span,
        int const& n_components) const
    {
        return ProcessLib::getIntegrationPointDataMaterialStateVariables(
            material_states_,
            &ConstitutiveRelations::MaterialStateData<
                DisplacementDim>::material_state_variables,
            get_values_span, n_components);
    }

    typename ConstitutiveRelations::SolidConstitutiveRelation<
        DisplacementDim>::MaterialStateVariables const&
    getMaterialStateVariablesAt(unsigned integration_point) const
    {
        return *material_states_[integration_point].material_state_variables;
    }

    static auto getReflectionDataForOutput()
    {
        using Self = LocalAssemblerInterface<DisplacementDim>;

        return ProcessLib::Reflection::reflectWithoutName(
            &Self::current_states_, &Self::output_data_);
    }

    TH2MProcessData<DisplacementDim>& process_data_;

    std::vector<typename ConstitutiveRelations::StatefulData<DisplacementDim>>
        current_states_;  // TODO maybe do not store but rather re-evaluate for
                          // state update
    std::vector<
        typename ConstitutiveRelations::StatefulDataPrev<DisplacementDim>>
        prev_states_;

    // Material state is special, because it contains both the current and the
    // old state.
    std::vector<ConstitutiveRelations::MaterialStateData<DisplacementDim>>
        material_states_;
    std::vector<ConstitutiveRelations::OutputData<DisplacementDim>>
        output_data_;

    NumLib::GenericIntegrationMethod const& integration_method_;
    MeshLib::Element const& element_;
    bool const is_axially_symmetric_;
    ConstitutiveRelations::SolidConstitutiveRelation<DisplacementDim> const&
        solid_material_;
};

}  // namespace TH2M
}  // namespace ProcessLib
