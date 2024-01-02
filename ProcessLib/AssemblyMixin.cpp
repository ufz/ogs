/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "AssemblyMixin.h"

#include <unordered_set>

#include "MeshLib/Utils/getOrCreateMeshProperty.h"
#include "NumLib/DOF/DOFTableUtil.h"

namespace
{
std::vector<std::reference_wrapper<MeshLib::PropertyVector<double>>>
createResiduumVectors(
    MeshLib::Mesh& mesh,
    std::vector<std::string> const& residuum_names,
    std::vector<std::reference_wrapper<ProcessLib::ProcessVariable>>
        pvs)
{
    auto const num_residua = residuum_names.size();

    std::vector<std::reference_wrapper<MeshLib::PropertyVector<double>>>
        residuum_vectors;
    residuum_vectors.reserve(num_residua);

    for (std::size_t i = 0; i < num_residua; ++i)
    {
        auto const& residuum_name = residuum_names[i];
        auto const& pv = pvs[i].get();

        residuum_vectors.emplace_back(*MeshLib::getOrCreateMeshProperty<double>(
            mesh, residuum_name, MeshLib::MeshItemType::Node,
            pv.getNumberOfGlobalComponents()));
    }

    return residuum_vectors;
}
}  // namespace

namespace ProcessLib
{
SubmeshAssemblyData::SubmeshAssemblyData(
    MeshLib::Mesh const& mesh,
    std::vector<std::reference_wrapper<MeshLib::PropertyVector<double>>>&&
        residuum_vectors)
    : bulk_element_ids{*MeshLib::bulkElementIDs(mesh)},
      bulk_node_ids{*MeshLib::bulkNodeIDs(mesh)},
      residuum_vectors{std::move(residuum_vectors)}
{
}

void AssemblyMixinBase::initializeAssemblyOnSubmeshes(
    const int process_id,
    MeshLib::Mesh& bulk_mesh,
    std::vector<std::reference_wrapper<MeshLib::Mesh>> const& submeshes,
    std::vector<std::string> const& residuum_names,
    std::vector<std::reference_wrapper<ProcessVariable>> const& pvs)
{
    DBUG("AssemblyMixinBase initializeSubmeshOutput().");

    auto const num_residua = residuum_names.size();

    if (pvs.size() != num_residua)
    {
        OGS_FATAL(
            "The number of passed residuum names does not match the number "
            "of process variables for process id {}: {} != {}",
            process_id, num_residua, pvs.size());
    }

    submesh_assembly_data_.reserve(submeshes.size());
    for (auto& mesh_ref : submeshes)
    {
        auto& mesh = mesh_ref.get();

        submesh_assembly_data_.emplace_back(
            mesh, createResiduumVectors(mesh, residuum_names, pvs));
    }

    residuum_vectors_bulk_ =
        createResiduumVectors(bulk_mesh, residuum_names, pvs);
}

void AssemblyMixinBase::updateActiveElements(
    ProcessLib::ProcessVariable const& pv)
{
    DBUG("AssemblyMixinBase updateActiveElements().");

    if (ids_state_ == ActiveElementIDsState::UNINITIALIZED)
    {
        updateActiveElementsImpl(pv);
        return;
    }

    ActiveElementIDsState const new_state =
        pv.getActiveElementIDs().empty()
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
    updateActiveElementsImpl(pv);
}

void AssemblyMixinBase::updateActiveElementsImpl(
    ProcessLib::ProcessVariable const& pv)
{
    DBUG("AssemblyMixinBase updateActiveElementsImpl().");

    auto const& active_element_ids = pv.getActiveElementIDs();

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
            // note: this copies the ID vector!
            sad.active_element_ids = sad.bulk_element_ids;
        }
    }
    else
    {
        // assemble each submesh on the intersection of the "global" active
        // elements and the submesh
        std::unordered_set<std::size_t> active_element_ids_set(
            active_element_ids.begin(), active_element_ids.end());

        for (auto& sad : submesh_assembly_data_)
        {
            auto& aeis = sad.active_element_ids;
            auto& beis = sad.bulk_element_ids;
            aeis.clear();
            aeis.reserve(beis.getNumberOfTuples());

            for (auto bei : beis)
            {
                if (active_element_ids_set.contains(bei))
                {
                    aeis.push_back(bei);
                }
            }
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
    GlobalVector const& rhs,
    NumLib::LocalToGlobalIndexMap const& local_to_global_index_map,
    SubmeshAssemblyData const& sad)
{
    for (std::size_t variable_id = 0; variable_id < sad.residuum_vectors.size();
         ++variable_id)
    {
        transformVariableFromGlobalVector(
            rhs, variable_id, local_to_global_index_map,
            sad.residuum_vectors[variable_id].get(), sad.bulk_node_ids,
            std::negate<double>());
    }
}

}  // namespace ProcessLib
