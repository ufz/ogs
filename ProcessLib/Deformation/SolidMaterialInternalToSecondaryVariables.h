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

template <typename SolidMaterial>
std::map<std::string,
         std::vector<std::pair<int, typename SolidMaterial::InternalVariable>>>
collectInternalVariables(
    std::map<int, std::unique_ptr<SolidMaterial>> const& solid_materials)
{
    assert(!solid_materials.empty());

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

    return internal_variables_by_name;
}

template <typename LocalAssemblerInterface, typename InternalVariable>
InternalVariable const* findInternalVariable(
    std::vector<std::pair<int, InternalVariable>> const& mat_iv_collection,
    bool const material_id_independent,
    LocalAssemblerInterface const& loc_asm)
{
    int const material_id =
        material_id_independent ? 0 : loc_asm.getMaterialID();

    auto const mat_iv_it = std::find_if(
        begin(mat_iv_collection), end(mat_iv_collection),
        [material_id](auto const& x) { return x.first == material_id; });
    if (mat_iv_it == end(mat_iv_collection))
    {
        // If local assembler does not provide correct solid material
        // model return empty vector, which will be ignored by the
        // extrapolation algorithm.
        return nullptr;
    }

    return &mat_iv_it->second;
}

template <typename SolidMaterial>
void forEachSolidMaterialInternalVariable(
    std::map<int, std::unique_ptr<SolidMaterial>> const& solid_materials,
    auto const& function)
{
    auto const internal_variables_by_name =
        collectInternalVariables(solid_materials);

    // Multiple material ids could be present but only one material for the
    // whole domain. In this case the choice of callbacks is independent of
    // local assembler's material id, and the material id is 0.
    // \see selectSolidConstitutiveRelation() for material id logic.
    bool const material_id_independent = solid_materials.size() == 1;

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

        function(name, num_components, mat_iv_collection,
                 material_id_independent);
    }
}

template <typename LocalAssemblerInterface, typename InternalVariable>
auto createCallback(
    std::vector<std::pair<int, InternalVariable>> const& mat_iv_collection,
    int const num_components,
    bool const material_id_independent)
{
    return [mat_iv_collection, num_components, material_id_independent](
               LocalAssemblerInterface const& loc_asm,
               const double /*t*/,
               std::vector<GlobalVector*> const& /*x*/,
               std::vector<
                   NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
               std::vector<double>& cache) -> std::vector<double> const&
    {
        cache.clear();

        auto const iv = findInternalVariable(mat_iv_collection,
                                             material_id_independent, loc_asm);

        if (iv == nullptr)
        {
            // If the material model for the present material group does not
            // have the requested internal variable, return an empty vector,
            // which will be ignored by the extrapolation algorithm.
            return cache;
        }

        auto const& fct = iv->getter;
        assert(num_components == iv->num_components);

        const unsigned num_int_pts = loc_asm.getNumberOfIntegrationPoints();
        assert(num_int_pts > 0);

        auto cache_mat = MathLib::createZeroedMatrix<Eigen::Matrix<
            double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
            cache, num_components, num_int_pts);

        // TODO avoid the heap allocation (one per finite element)
        std::vector<double> cache_column(num_int_pts);

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
}

template <typename LocalAssemblerInterface, typename SolidMaterial,
          typename AddSecondaryVariableCallback>
void solidMaterialInternalToSecondaryVariables(
    std::map<int, std::unique_ptr<SolidMaterial>> const& solid_materials,
    AddSecondaryVariableCallback const& add_secondary_variable)
{
    auto register_secondary_variable =
        [&add_secondary_variable](
            std::string name, char const num_components,
            std::vector<
                std::pair<int, typename SolidMaterial::InternalVariable>> const&
                mat_iv_collection,
            bool const material_id_independent)
    {
        DBUG("Registering internal variable {:s}.", name);

        add_secondary_variable(
            name, num_components,
            createCallback<LocalAssemblerInterface>(
                mat_iv_collection, num_components, material_id_independent));
    };

    forEachSolidMaterialInternalVariable(solid_materials,
                                         register_secondary_variable);
}

template <typename LocalAssemblerInterface, typename InternalVariable>
auto createCallbackForIpWriter(
    std::vector<std::pair<int, InternalVariable>> const& mat_iv_collection,
    int const num_components,
    bool const material_id_independent)
{
    return [mat_iv_collection, num_components, material_id_independent](
               LocalAssemblerInterface const& loc_asm) -> std::vector<double>
    {
        auto const iv = findInternalVariable(mat_iv_collection,
                                             material_id_independent, loc_asm);
        if (iv == nullptr)
        {
            // If the material model for the present material group does not
            // have the requested internal variable, return an empty vector,
            // which will be ignored by the integration point writer.
            return {};
        }

        auto const& fct = iv->reference;
        assert(num_components == iv->num_components);

        return loc_asm.getMaterialStateVariableInternalState(fct,
                                                             num_components);
    };
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
    auto add_integration_point_writer =
        [&local_assemblers, &integration_point_writer, integration_order](
            std::string name, char const num_components,
            std::vector<
                std::pair<int, typename SolidMaterial::InternalVariable>> const&
                mat_iv_collection,
            bool const material_id_independent)
    {
        DBUG("Creating integration point writer for  internal variable {:s}.",
             name);

        integration_point_writer.emplace_back(
            std::make_unique<MeshLib::IntegrationPointWriter>(
                "material_state_variable_" + name + "_ip", num_components,
                integration_order, local_assemblers,
                createCallbackForIpWriter<LocalAssemblerInterface>(
                    mat_iv_collection, num_components,
                    material_id_independent)));
    };

    forEachSolidMaterialInternalVariable(solid_materials,
                                         add_integration_point_writer);
}
}  // namespace ProcessLib::Deformation
