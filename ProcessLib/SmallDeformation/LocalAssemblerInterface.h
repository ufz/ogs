/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>

#include "ConstitutiveRelations/MaterialState.h"
#include "MaterialLib/SolidModels/MechanicsBase.h"
#include "MaterialLib/SolidModels/SelectSolidConstitutiveRelation.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/Integration/GenericIntegrationMethod.h"
#include "ProcessLib/Deformation/MaterialForces.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "ProcessLib/Utils/SetOrGetIntegrationPointData.h"
#include "SmallDeformationProcessData.h"

namespace ProcessLib
{
namespace SmallDeformation
{
template <int DisplacementDim>
struct SmallDeformationLocalAssemblerInterface
    : public ProcessLib::LocalAssemblerInterface,
      public MaterialForcesInterface,
      public NumLib::ExtrapolatableElement
{
    SmallDeformationLocalAssemblerInterface(
        MeshLib::Element const& e,
        NumLib::GenericIntegrationMethod const& integration_method,
        bool const is_axially_symmetric,
        SmallDeformationProcessData<DisplacementDim>& process_data)
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

        material_states_.reserve(n_integration_points);
        for (unsigned ip = 0; ip < n_integration_points; ++ip)
        {
            material_states_.emplace_back(
                solid_material_.createMaterialStateVariables());
        }
    }
    virtual std::size_t setIPDataInitialConditions(
        std::string const& name, double const* values,
        int const integration_order) = 0;

    virtual std::vector<double> const& getIntPtFreeEnergyDensity(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> getSigma() const = 0;
    virtual std::vector<double> const& getIntPtSigma(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtEpsilon(
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
            typename MaterialLib::Solids::MechanicsBase<DisplacementDim>::
                MaterialStateVariables&)> const& get_values_span,
        int const& n_components) const
    {
        return ProcessLib::getIntegrationPointDataMaterialStateVariables(
            material_states_,
            &MaterialStateData<DisplacementDim>::material_state_variables,
            get_values_span, n_components);
    }

    typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::MaterialStateVariables const&
    getMaterialStateVariablesAt(unsigned integration_point) const
    {
        return *material_states_[integration_point].material_state_variables;
    }

protected:
    SmallDeformationProcessData<DisplacementDim>& process_data_;

    // Material state is special, because it contains both the current and the
    // old state.
    std::vector<MaterialStateData<DisplacementDim>> material_states_;

    NumLib::GenericIntegrationMethod const& integration_method_;
    MeshLib::Element const& element_;
    bool const is_axially_symmetric_;
    MaterialLib::Solids::MechanicsBase<DisplacementDim> const& solid_material_;
};

}  // namespace SmallDeformation
}  // namespace ProcessLib
