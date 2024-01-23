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

#include <Eigen/Eigenvalues>
#include <vector>

#include "ConstitutiveRelations/ConstitutiveData.h"
#include "ConstitutiveRelations/ConstitutiveSetting.h"
#include "ConstitutiveRelations/MaterialState.h"
#include "MaterialForces.h"
#include "MaterialLib/SolidModels/MechanicsBase.h"
#include "MaterialLib/SolidModels/SelectSolidConstitutiveRelation.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/Integration/GenericIntegrationMethod.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "ProcessLib/Reflection/ReflectionSetIPData.h"
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
            this->current_states_[ip].stress_data.sigma.noalias() =
                MathLib::KelvinVector::KelvinVectorType<
                    DisplacementDim>::Zero();
        }
    }
    /// Returns number of read integration points.
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

            WARN(
                "Could not find variable {:s} in solid material model's "
                "internal variables.",
                name);
            return 0;
        }

        // TODO this logic could be pulled out of the local assembler into the
        // process. That might lead to a slightly better performance due to less
        // string comparisons.
        return ProcessLib::Reflection::reflectSetIPData<DisplacementDim>(
            name, values, current_states_);
    }

    void computeSecondaryVariableConcrete(
        double const /*t*/, double const /*dt*/, Eigen::VectorXd const& /*x*/,
        Eigen::VectorXd const& /*x_prev*/) override
    {
        int const elem_id = element_.getID();
        ParameterLib::SpatialPosition x_position;
        x_position.setElementID(elem_id);
        unsigned const n_integration_points =
            integration_method_.getNumberOfPoints();

        auto sigma_sum = MathLib::KelvinVector::tensorToKelvin<DisplacementDim>(
            Eigen::Matrix<double, 3, 3>::Zero());

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            x_position.setIntegrationPoint(ip);
            auto const& sigma = current_states_[ip].stress_data.sigma;
            sigma_sum += sigma;
        }

        Eigen::Matrix<double, 3, 3, 0, 3, 3> const sigma_avg =
            MathLib::KelvinVector::kelvinVectorToTensor(sigma_sum) /
            n_integration_points;

        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 3, 3>> e_s(
            sigma_avg);

        Eigen::Map<Eigen::Vector3d>(
            &(*process_data_.principal_stress_values)[elem_id * 3], 3) =
            e_s.eigenvalues();

        auto eigen_vectors = e_s.eigenvectors();

        for (auto i = 0; i < 3; i++)
        {
            Eigen::Map<Eigen::Vector3d>(
                &(*process_data_.principal_stress_vector[i])[elem_id * 3], 3) =
                eigen_vectors.col(i);
        }
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

    static auto getReflectionDataForOutput()
    {
        using Self = SmallDeformationLocalAssemblerInterface<DisplacementDim>;

        return ProcessLib::Reflection::reflectWithoutName(
            &Self::current_states_, &Self::output_data_);
    }

protected:
    SmallDeformationProcessData<DisplacementDim>& process_data_;

    // Material state is special, because it contains both the current and the
    // old state.
    std::vector<MaterialStateData<DisplacementDim>> material_states_;
    std::vector<typename ConstitutiveRelations::StatefulData<DisplacementDim>>
        current_states_;  // TODO maybe do not store but rather re-evaluate for
                          // state update
    std::vector<
        typename ConstitutiveRelations::StatefulDataPrev<DisplacementDim>>
        prev_states_;
    std::vector<typename ConstitutiveRelations::OutputData<DisplacementDim>>
        output_data_;

    NumLib::GenericIntegrationMethod const& integration_method_;
    MeshLib::Element const& element_;
    bool const is_axially_symmetric_;
    MaterialLib::Solids::MechanicsBase<DisplacementDim> const& solid_material_;
    typename ConstitutiveRelations::ConstitutiveSetting<DisplacementDim>
        constitutive_setting;
};

}  // namespace SmallDeformation
}  // namespace ProcessLib
