// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "AssemblyMixin.h"

#include <range/v3/algorithm/sort.hpp>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/for_each.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/zip.hpp>

#include "MeshLib/Utils/getOrCreateMeshProperty.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "Process.h"

namespace
{

void checkResiduumNamesVsProcessVariables(
    std::vector<std::vector<std::string>> const& per_process_residuum_names,
    std::vector<
        std::vector<std::reference_wrapper<ProcessLib::ProcessVariable>>> const&
        per_process_pvs)
{
    if (per_process_pvs.size() != per_process_residuum_names.size())
    {
        OGS_FATAL(
            "The number of passed residuum names ({}) does not match the "
            "number of processes ({}).",
            per_process_residuum_names.size(), per_process_pvs.size());
    }

    auto check_sizes = [](std::size_t const process_id, auto const& rns_pvs)
    {
        auto const& [rns, pvs] = rns_pvs;

        if (rns.size() != pvs.size())
        {
            OGS_FATAL(
                "The number of passed residuum names ({}) does not match the "
                "number of process variables ({}) for process {}.",
                rns.size(), pvs.size(), process_id);
        }
    };
    for (auto const& [process_id, rns_pvs] :
         ranges::views::zip(per_process_residuum_names, per_process_pvs) |
             ranges::views::enumerate)
    {
        check_sizes(process_id, rns_pvs);
    }
}

std::vector<
    std::vector<std::reference_wrapper<MeshLib::PropertyVector<double>>>>
createResiduumVectors(
    MeshLib::Mesh& mesh,
    std::vector<std::vector<std::string>> const& per_process_residuum_names,
    std::vector<
        std::vector<std::reference_wrapper<ProcessLib::ProcessVariable>>> const&
        per_process_pvs)
{
    checkResiduumNamesVsProcessVariables(per_process_residuum_names,
                                         per_process_pvs);

    auto create_mesh_property_for_residuum =
        [&mesh](std::pair<std::string const&,
                          ProcessLib::ProcessVariable const&> const&
                    residuum_name_process_variable)
        -> std::reference_wrapper<MeshLib::PropertyVector<double>>
    {
        auto const& [rn, pv] = residuum_name_process_variable;
        return *MeshLib::getOrCreateMeshProperty<double>(
            mesh, rn, MeshLib::MeshItemType::Node,
            pv.getNumberOfGlobalComponents());
    };

    auto create_mesh_properties =
        [&create_mesh_property_for_residuum](auto const& rns_pvs)
    {
        auto const& [rns, pvs] = rns_pvs;
        return ranges::views::zip(rns, pvs) |
               ranges::views::transform(create_mesh_property_for_residuum) |
               ranges::to<std::vector<
                   std::reference_wrapper<MeshLib::PropertyVector<double>>>>;
    };

    return ranges::views::zip(per_process_residuum_names, per_process_pvs) |
           ranges::views::transform(create_mesh_properties) |
           ranges::to<std::vector<std::vector<
               std::reference_wrapper<MeshLib::PropertyVector<double>>>>>;
}
}  // namespace

namespace ProcessLib
{

AssemblyMixinBase::AssemblyMixinBase(
    AbstractJacobianAssembler& jacobian_assembler,
    bool const is_linear,
    bool const use_monolithic_scheme)
    : pvma_{jacobian_assembler}, is_linear_{is_linear}
{
    if (is_linear && !use_monolithic_scheme)
    {
        OGS_FATAL(
            "You requested to assemble only once in combination with staggered "
            "coupling. This use case is not yet implemented.");
    }

    if (is_linear_)
    {
        WARN(
            "You specified that the process simulated by OGS is linear. With "
            "this optimization, the system is assembled only once, and the "
            "non-linear solver performs a single iteration per time step. No "
            "non-linearities will be resolved, and OGS will not detect whether "
            "any exist. It is your responsibility to ensure that the assembled "
            "systems are indeed linear; there is no safety net.\n"
            "\n"
            "Be aware that non-linearities may arise in subtle ways. For "
            "example, in an HT problem with a steady-state velocity field, if "
            "the initial pressure field does not match the steady-state "
            "condition (e.g., initialized to zero), then the steady-state "
            "velocity field only emerges during the first iteration—making the "
            "system effectively non-linear under this optimization.");
    }
}

void AssemblyMixinBase::initializeAssemblyOnSubmeshes(
    MeshLib::Mesh& bulk_mesh,
    std::vector<std::reference_wrapper<MeshLib::Mesh>> const& submeshes,
    std::vector<std::vector<std::string>> const& residuum_names,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>> const&
        pvs)
{
    DBUG("AssemblyMixinBase initializeSubmeshOutput().");

    auto const num_submeshes = submeshes.size();

    if (num_submeshes == 0)
    {
        bulk_mesh_assembly_data_.emplace(
            createResiduumVectors(bulk_mesh, residuum_names, pvs));
        return;
    }

    submesh_assembly_data_.reserve(num_submeshes);
    for (auto& submesh_ref : submeshes)
    {
        auto& submesh = submesh_ref.get();

        submesh_assembly_data_.emplace_back(
            submesh, createResiduumVectors(submesh, residuum_names, pvs));
    }

    submesh_matrix_cache_.resize(num_submeshes);
}

void AssemblyMixinBase::updateActiveElements(Process const& process)
{
    DBUG("AssemblyMixin updateActiveElements().");

    if (ids_state_ == ActiveElementIDsState::UNINITIALIZED)
    {
        updateActiveElementsImpl(process);
        return;
    }

    ActiveElementIDsState const new_state =
        process.getActiveElementIDs().empty()
            ? ActiveElementIDsState::NO_DEACTIVATED_SUBDOMAINS
            : ActiveElementIDsState::HAS_DEACTIVATED_SUBDOMAINS;

    if (ids_state_ == ActiveElementIDsState::NO_DEACTIVATED_SUBDOMAINS &&
        new_state == ActiveElementIDsState::NO_DEACTIVATED_SUBDOMAINS)
    {
        // no update necessary
        return;
    }

    // updating the active elements is necessary in all remaining cases, because
    // of the following transitions between old and new states (deactivated
    // subdomains present? yes/no; old state -> new state):
    // * no -> yes - now there are deactivated subdomains
    // * yes -> no - no deactivated subdomains anymore
    // * yes -> yes - deactivated subdomains might have changed
    updateActiveElementsImpl(process);
}

void AssemblyMixinBase::updateActiveElementsImpl(Process const& process)
{
    DBUG("AssemblyMixinBase updateActiveElementsImpl().");

    auto const& active_element_ids = process.getActiveElementIDs();

    ActiveElementIDsState const new_state =
        active_element_ids.empty()
            ? ActiveElementIDsState::NO_DEACTIVATED_SUBDOMAINS
            : ActiveElementIDsState::HAS_DEACTIVATED_SUBDOMAINS;

    if (new_state == ActiveElementIDsState::NO_DEACTIVATED_SUBDOMAINS)
    {
        // no deactivated subdomains, assemble on all elements as specified
        // in the bulk_element_ids of the submeshes
        for (auto& sad : submesh_assembly_data_)
        {
            sad.setAllElementsActive();
        }
        if (bulk_mesh_assembly_data_)
        {
            bulk_mesh_assembly_data_->setAllElementsActive();
        }
    }
    else
    {
        auto active_element_ids_copy{active_element_ids};
        // assemble each submesh on the intersection of the "global" active
        // elements and the submesh
        ranges::sort(active_element_ids_copy);

        for (auto& sad : submesh_assembly_data_)
        {
            sad.setElementSelectionActive(active_element_ids_copy);
        }
        if (bulk_mesh_assembly_data_)
        {
            bulk_mesh_assembly_data_->setElementSelectionActive(
                active_element_ids_copy);
        }
    }

    ids_state_ = new_state;
}

void AssemblyMixinBase::copyResiduumVectorsToBulkMesh(
    GlobalVector const& rhs,
    NumLib::LocalToGlobalIndexMap const& local_to_global_index_map,
    std::vector<std::reference_wrapper<MeshLib::PropertyVector<double>>>
        residuum_vectors)
{
    for (std::size_t variable_id = 0; variable_id < residuum_vectors.size();
         ++variable_id)
    {
        transformVariableFromGlobalVector(
            rhs, variable_id, local_to_global_index_map,
            residuum_vectors[variable_id].get(), std::negate<double>());
    }
}

void AssemblyMixinBase::copyResiduumVectorsToSubmesh(
    int const process_id,
    GlobalVector const& rhs,
    NumLib::LocalToGlobalIndexMap const& local_to_global_index_map,
    ProcessLib::Assembly::SubmeshAssemblyData const& sad)
{
    auto const& residuum_vectors = sad.residuum_vectors[process_id];
    for (std::size_t variable_id = 0; variable_id < residuum_vectors.size();
         ++variable_id)
    {
        transformVariableFromGlobalVector(
            rhs, variable_id, local_to_global_index_map,
            residuum_vectors[variable_id].get(),
            std::span<std::size_t const>{sad.bulk_node_ids},
            std::negate<double>());
    }
}

}  // namespace ProcessLib
