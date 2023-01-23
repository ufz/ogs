/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <map>
#include <memory>

#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "MathLib/LinAlg/GlobalMatrixVectorTypes.h"
#include "MeshLib/Utils/IntegrationPointWriter.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"

namespace ProcessLib::Deformation
{
template <typename LocalAssemblerInterface, typename SolidMaterial,
          typename AddSecondaryVariableCallback>
void solidMaterialInternalToSecondaryVariables(
    std::map<int, std::unique_ptr<SolidMaterial>> const& solid_materials,
    AddSecondaryVariableCallback const& add_secondary_variable)
{
    // For each name of an internal variable collect all solid material/
    // internal variable pairs.
    std::map<
        std::string,
        std::vector<std::pair<int, typename SolidMaterial::InternalVariable>>>
        internal_variables_by_name;
    for (auto const& [material_id, solid_material] : solid_materials)
    {
        auto const& internal_variables = solid_material->getInternalVariables();
        for (auto const& iv : internal_variables)
        {
            internal_variables_by_name[iv.name].push_back({material_id, iv});
        }
    }

    // Create *single* callback passing all solid materials to it. Choose
    // correct solid material based on the local assembler's solid material in
    // the callback.
    for (auto const& [name, mat_iv_collection] : internal_variables_by_name)
    {
        assert(!mat_iv_collection.empty());
        auto const num_components =
            mat_iv_collection.front().second.num_components;

        // Check that the number of components is equal for all materials.
        if (!std::all_of(
                begin(mat_iv_collection), end(mat_iv_collection),
                [num_components](auto const& mat_iv)
                { return mat_iv.second.num_components == num_components; }))
        {
            OGS_FATAL(
                "Not for all material ids the secondary variable '{:s}' has "
                "{:d} components.",
                name, num_components);
        }

        DBUG("Registering internal variable {:s}.", name);

        auto callback =
            [mat_iv_collection = mat_iv_collection, num_components](
                LocalAssemblerInterface const& loc_asm,
                const double /*t*/,
                std::vector<GlobalVector*> const& /*x*/,
                std::vector<
                    NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
                std::vector<double>& cache) -> std::vector<double> const&
        {
            cache.clear();

            const unsigned num_int_pts = loc_asm.getNumberOfIntegrationPoints();
            assert(num_int_pts > 0);

            int const material_id = loc_asm.getMaterialID();

            auto const mat_iv =
                std::find_if(begin(mat_iv_collection), end(mat_iv_collection),
                             [material_id](auto const& x)
                             { return x.first == material_id; });
            if (mat_iv == end(mat_iv_collection))
            {
                // If local assembler does not provide correct solid material
                // model return empty vector, which will be ignored by the
                // extrapolation algorithm.
                return cache;
            }

            auto cache_mat = MathLib::createZeroedMatrix<Eigen::Matrix<
                double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
                cache, num_components, num_int_pts);

            // TODO avoid the heap allocation (one per finite element)
            std::vector<double> cache_column(num_int_pts);

            auto const& iv = mat_iv->second;
            auto const& fct = iv.getter;

            assert(material_id == mat_iv->first);
            assert(num_components == iv.num_components);

            for (unsigned i = 0; i < num_int_pts; ++i)
            {
                auto const& state = loc_asm.getMaterialStateVariablesAt(i);

                auto const& int_pt_values = fct(state, cache_column);
                assert(int_pt_values.size() ==
                       static_cast<std::size_t>(num_components));
                auto const int_pt_values_vec = MathLib::toVector(int_pt_values);

                cache_mat.col(i).noalias() = int_pt_values_vec;
            }

            return cache;
        };

        add_secondary_variable(name, num_components, std::move(callback));
    }
}

template <typename LocalAssemblerInterface, typename SolidMaterial>
void solidMaterialInternalVariablesToIntegrationPointWriter(
    std::map<int, std::unique_ptr<SolidMaterial>> const& solid_materials,
    std::vector<std::unique_ptr<LocalAssemblerInterface>> const&
        local_assemblers,
    std::vector<std::unique_ptr<MeshLib::IntegrationPointWriter>>&
        integration_point_writer,
    int const integration_order)
{
    // Collect the internal variables for all solid materials.
    std::vector<typename SolidMaterial::InternalVariable> internal_variables;
    for (auto const& solid_material : solid_materials)
    {
        auto const ivs = solid_material.second->getInternalVariables();
        copy(begin(ivs), end(ivs), back_inserter(internal_variables));
    }

    // Create integration point writers for each of the internal variables.
    for (auto const& iv : internal_variables)
    {
        DBUG("Creating integration point writer for  internal variable {:s}.",
             iv.name);

        integration_point_writer.emplace_back(
            std::make_unique<MeshLib::IntegrationPointWriter>(
                "material_state_variable_" + iv.name + "_ip", iv.num_components,
                integration_order, local_assemblers,
                &LocalAssemblerInterface::getMaterialStateVariableInternalState,
                iv.reference, iv.num_components));
    }
}
}  // namespace ProcessLib::Deformation
