/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <map>
#include <memory>

#include "MathLib/LinAlg/Eigen/EigenMapTools.h"

namespace ProcessLib::Deformation
{
template <typename LocalAssemblerInterface,
          typename AddSecondaryVariableCallback, int DisplacementDim>
void solidMaterialInternalToSecondaryVariables(
    std::map<int, std::unique_ptr<MaterialLib::Solids::MechanicsBase<
                      DisplacementDim>>> const& solid_materials,
    AddSecondaryVariableCallback const& add_secondary_variable)
{
    // Collect the internal variables for all solid materials.
    std::vector<typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::InternalVariable>
        internal_variables;
    for (auto const& material_id__solid_material : solid_materials)
    {
        auto const variables =
            material_id__solid_material.second->getInternalVariables();
        copy(begin(variables), end(variables),
             back_inserter(internal_variables));
    }

    // Register the internal variables.
    for (auto const& internal_variable : internal_variables)
    {
        auto const& name = internal_variable.name;
        auto const& fct = internal_variable.getter;
        auto const num_components = internal_variable.num_components;
        DBUG("Registering internal variable {:s}.", name);

        auto getIntPtValues =
            [fct, num_components](
                LocalAssemblerInterface const& loc_asm,
                const double /*t*/,
                std::vector<GlobalVector*> const& /*x*/,
                std::vector<
                    NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
                std::vector<double>& cache) -> std::vector<double> const& {
            const unsigned num_int_pts = loc_asm.getNumberOfIntegrationPoints();

            cache.clear();
            auto cache_mat = MathLib::createZeroedMatrix<Eigen::Matrix<
                double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
                cache, num_components, num_int_pts);

            // TODO avoid the heap allocation (one per finite element)
            std::vector<double> cache_column(num_int_pts);

            for (unsigned i = 0; i < num_int_pts; ++i)
            {
                auto const& state = loc_asm.getMaterialStateVariablesAt(i);

                auto const& int_pt_values = fct(state, cache_column);
                assert(int_pt_values.size() == num_components);
                auto const int_pt_values_vec = MathLib::toVector(int_pt_values);

                cache_mat.col(i).noalias() = int_pt_values_vec;
            }

            return cache;
        };

        add_secondary_variable(name, num_components, std::move(getIntPtValues));
    }
}

template <typename LocalAssemblerInterface, typename AddIntegrationPointWriter,
          int DisplacementDim>
void solidMaterialInternalVariablesToIntegrationPointWriter(
    std::map<int, std::unique_ptr<MaterialLib::Solids::MechanicsBase<
                      DisplacementDim>>> const& solid_materials,
    std::vector<std::unique_ptr<LocalAssemblerInterface>> const&
        local_assemblers,
    AddIntegrationPointWriter const& add_integration_point_writer)
{
    // Collect the internal variables for all solid materials.
    std::vector<typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::InternalVariable>
        internal_variables;
    for (auto const& solid_material : solid_materials)
    {
        auto const ivs = solid_material.second->getInternalVariables();
        copy(begin(ivs), end(ivs), back_inserter(internal_variables));
    }

    // Create integration point writers for each of the internal variables.
    for (auto const& internal_variable : internal_variables)
    {
        auto const& name = internal_variable.name;
        auto const& fct = internal_variable.reference;
        auto const n_components = internal_variable.num_components;
        DBUG("Creating integration point writer for  internal variable {:s}.",
             name);

        auto get_integration_point_data = [fct, n_components,
                                           &local_assemblers]() {
            // Result containing integration point data for each local
            // assembler.
            std::vector<std::vector<double>> result;
            result.reserve(local_assemblers.size());

            std::transform(begin(local_assemblers), end(local_assemblers),
                           back_inserter(result),
                           [fct, n_components](auto const& local_assembler) {
                               return local_assembler
                                   ->getMaterialStateVariableInternalState(
                                       fct, n_components);
                           });

            return result;
        };

        add_integration_point_writer(name, n_components,
                                     get_integration_point_data);
    }
}
}  // namespace ProcessLib::Deformation
