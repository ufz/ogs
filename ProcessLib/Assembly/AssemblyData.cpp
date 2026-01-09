// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "AssemblyData.h"

#include <range/v3/algorithm/is_sorted.hpp>
#include <range/v3/algorithm/set_algorithm.hpp>
#include <range/v3/algorithm/sort.hpp>

namespace ProcessLib::Assembly
{
std::shared_ptr<std::vector<std::size_t> const>
CommonAssemblyData::activeElementIDsSorted(
    std::vector<std::size_t> const* const sorted_element_subset) const
{
    if (!sorted_element_subset)
    {
        return sorted_active_element_ids_;
    }

    assert(ranges::is_sorted(*sorted_element_subset));

    if (areAllElementsActive())
    {
        // We do not own sorted_element_subset, so the returned data must not be
        // deleted by us.
        auto do_not_delete = [](auto*) {};
        return {sorted_element_subset, do_not_delete};
    }

    // Final case: there are both some deactivated elements and a further
    // subselection by sorted_element_subset.

    // TODO frequent re-computations might occur here. Optimization potential.
    auto aeis = std::make_shared<std::vector<std::size_t>>();
    aeis->reserve(std::min(sorted_active_element_ids_->size(),
                           sorted_element_subset->size()));

    ranges::set_intersection(*sorted_active_element_ids_,
                             *sorted_element_subset, std::back_inserter(*aeis));
    return aeis;
}

void BulkMeshAssemblyData::setAllElementsActive()
{
    // A nullptr means all elements are active.
    sorted_active_element_ids_.reset();
}

void BulkMeshAssemblyData::setElementSelectionActive(
    std::vector<std::size_t> const& sorted_active_element_ids_whole_mesh)
{
    assert(ranges::is_sorted(sorted_active_element_ids_whole_mesh));

    sorted_active_element_ids_ =
        std::make_shared<std::vector<std::size_t> const>(
            sorted_active_element_ids_whole_mesh);
}

SubmeshAssemblyData::SubmeshAssemblyData(
    MeshLib::Mesh const& submesh,
    std::vector<
        std::vector<std::reference_wrapper<MeshLib::PropertyVector<double>>>>&&
        residuum_vectors)
    : CommonAssemblyData{std::move(residuum_vectors)},
      bulk_node_ids{*MeshLib::bulkNodeIDs(submesh)},
      bulk_element_ids{*MeshLib::bulkElementIDs(submesh)}
{
}

void SubmeshAssemblyData::setAllElementsActive()
{
    std::vector<std::size_t> aeis;
    aeis.insert(aeis.end(), bulk_element_ids.begin(), bulk_element_ids.end());
    ranges::sort(aeis);
    sorted_active_element_ids_ =
        std::make_shared<std::vector<std::size_t> const>(std::move(aeis));
}

void SubmeshAssemblyData::setElementSelectionActive(
    std::vector<std::size_t> const& sorted_active_element_ids_whole_mesh)
{
    assert(ranges::is_sorted(sorted_active_element_ids_whole_mesh));

    std::vector<std::size_t> aeis;

    ranges::set_intersection(bulk_element_ids,
                             sorted_active_element_ids_whole_mesh,
                             std::back_inserter(aeis));

    sorted_active_element_ids_ =
        std::make_shared<std::vector<std::size_t> const>(std::move(aeis));
}
}  // namespace ProcessLib::Assembly
