/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
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

    // Check for each internal variable name, that the number of components is
    // equal for all materials.
    for (auto const& [name, mat_iv_collection] : internal_variables_by_name)
    {
        if (mat_iv_collection.empty())
        {
            continue;
        }
        auto const num_components =
            mat_iv_collection.front().second.num_components;

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
    }

    return internal_variables_by_name;
}

template <typename InternalVariable>
struct InternalVariablesCollection
{
    template <typename LocalAssemblerInterface>
    InternalVariable const* findInternalVariable(
        LocalAssemblerInterface const& loc_asm) const
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

        assert(num_components == mat_iv_it->second.num_components);

        return &mat_iv_it->second;
    }

    std::vector<std::pair<int, InternalVariable>> mat_iv_collection;
    int num_components;
    bool material_id_independent;
};

template <typename Pair>
InternalVariablesCollection(std::vector<Pair>&&, int, bool)
    -> InternalVariablesCollection<typename Pair::second_type>;

template <typename SolidMaterial>
void forEachSolidMaterialInternalVariable(
    std::map<int, std::unique_ptr<SolidMaterial>> const& solid_materials,
    auto const& function)
{
    auto internal_variables_by_name = collectInternalVariables(solid_materials);

    // Multiple material ids could be present but only one material for the
    // whole domain. In this case the choice of callbacks is independent of
    // local assembler's material id, and the material id is 0.
    // \see selectSolidConstitutiveRelation() for material id logic.
    bool const material_id_independent = solid_materials.size() == 1;

    // Create *single* callback passing all solid materials to it. Choose
    // correct solid material based on the local assembler's solid material in
    // the callback.
    for (auto&& [name, mat_iv_collection] : internal_variables_by_name)
    {
        auto const num_components =
            mat_iv_collection.front().second.num_components;

        function(name,
                 InternalVariablesCollection{std::move(mat_iv_collection),
                                             num_components,
                                             material_id_independent});
    }
}

template <typename InternalVariable>
class CollectIntegrationPointDataForExtrapolation final
{
public:
    explicit CollectIntegrationPointDataForExtrapolation(
        InternalVariablesCollection<InternalVariable>&&
            internal_variables_collection)
        : internal_variables_collection_(
              std::move(internal_variables_collection))
    {
    }

    template <typename LocalAssemblerInterface>
    std::vector<double> const& operator()(
        LocalAssemblerInterface const& loc_asm,
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const
    {
        cache.clear();

        auto const* const iv =
            internal_variables_collection_.findInternalVariable(loc_asm);

        if (iv == nullptr)
        {
            // If the material model for the present material group does not
            // have the requested internal variable, return an empty vector,
            // which will be ignored by the extrapolation algorithm.
            return cache;
        }

        auto const& fct = iv->getter;

        const unsigned num_int_pts = loc_asm.getNumberOfIntegrationPoints();
        assert(num_int_pts > 0);

        auto cache_mat = MathLib::createZeroedMatrix<Eigen::Matrix<
            double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
            cache, internal_variables_collection_.num_components, num_int_pts);

        // TODO avoid the heap allocation (one per finite element)
        std::vector<double> cache_column(num_int_pts);

        for (unsigned i = 0; i < num_int_pts; ++i)
        {
            auto const& state = loc_asm.getMaterialStateVariablesAt(i);

            auto const& int_pt_values = fct(state, cache_column);
            assert(int_pt_values.size() ==
                   static_cast<std::size_t>(
                       internal_variables_collection_.num_components));
            auto const int_pt_values_vec = MathLib::toVector(int_pt_values);

            cache_mat.col(i).noalias() = int_pt_values_vec;
        }

        return cache;
    }

private:
    InternalVariablesCollection<InternalVariable>
        internal_variables_collection_;
};

template <typename LocalAssemblerInterface, typename SolidMaterial,
          typename AddSecondaryVariableCallback>
void solidMaterialInternalToSecondaryVariables(
    std::map<int, std::unique_ptr<SolidMaterial>> const& solid_materials,
    AddSecondaryVariableCallback const& add_secondary_variable)
{
    auto register_secondary_variable =
        [&add_secondary_variable](
            std::string const& name,
            InternalVariablesCollection<
                typename SolidMaterial::InternalVariable>&&
                internal_variables_collection)
    {
        DBUG("Registering internal variable {:s}.", name);

        add_secondary_variable(name,
                               internal_variables_collection.num_components,
                               CollectIntegrationPointDataForExtrapolation{
                                   std::move(internal_variables_collection)});
    };

    forEachSolidMaterialInternalVariable(solid_materials,
                                         register_secondary_variable);
}

template <typename InternalVariable>
class CollectIntegrationPointDataForIpWriter final
{
public:
    explicit CollectIntegrationPointDataForIpWriter(
        InternalVariablesCollection<InternalVariable>&&
            internal_variables_collection)
        : internal_variables_collection_(internal_variables_collection)
    {
    }

    template <typename LocalAssemblerInterface>
    std::vector<double> operator()(LocalAssemblerInterface const& loc_asm) const
    {
        auto const iv =
            internal_variables_collection_.findInternalVariable(loc_asm);
        if (iv == nullptr)
        {
            // If the material model for the present material group does not
            // have the requested internal variable, return an empty vector,
            // which will be ignored by the integration point writer.
            return {};
        }

        auto const& fct = iv->reference;

        return loc_asm.getMaterialStateVariableInternalState(
            fct, internal_variables_collection_.num_components);
    };

private:
    InternalVariablesCollection<InternalVariable>
        internal_variables_collection_;
};

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
            std::string const& name,
            InternalVariablesCollection<
                typename SolidMaterial::InternalVariable>&&
                internal_variables_collection)
    {
        DBUG("Creating integration point writer for  internal variable {:s}.",
             name);

        integration_point_writer.emplace_back(
            std::make_unique<MeshLib::IntegrationPointWriter>(
                "material_state_variable_" + name + "_ip",
                internal_variables_collection.num_components, integration_order,
                local_assemblers,
                CollectIntegrationPointDataForIpWriter{
                    std::move(internal_variables_collection)}));
    };

    forEachSolidMaterialInternalVariable(solid_materials,
                                         add_integration_point_writer);
}
}  // namespace ProcessLib::Deformation
