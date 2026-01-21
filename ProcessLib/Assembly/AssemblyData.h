// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "MeshLib/Mesh.h"
#include "MeshLib/PropertyVector.h"

namespace ProcessLib::Assembly
{
//! Data necessary for global equation system assembly.
//!
//! Base of both bulk mesh and submesh assembly data.
struct CommonAssemblyData
{
    explicit CommonAssemblyData(
        std::vector<std::vector<
            std::reference_wrapper<MeshLib::PropertyVector<double>>>>&&
            residuum_vectors)
        : residuum_vectors{std::move(residuum_vectors)}
    {
    }

    /*! Returns active element ids.
     *
     * \pre
     * setAllElementsActive() or setElementSelectionActive() must have been
     * called before.
     *
     * If sorted_element_subset is not nullptr, only the intersection between
     * the previously set values and sorted_element_subset is returned.
     * sorted_element_subset must be sorted in ascending order in that case. At
     * the moment (Nov 2025) this is used for assembly optimizations in the
     * HeatTransportBHE process.
     *
     * \return
     * The returned range of element IDs is sorted in ascending order. (NOTE:
     * That's specified for std::ranges::set_intersection only, I hope it is
     * true for range-v3, too.)
     * If this method returns nullptr it means that all elements are active.
     */
    std::shared_ptr<std::vector<std::size_t> const> activeElementIDsSorted(
        std::vector<std::size_t> const* const sorted_element_subset) const;

    //! Residuum vectors for each process ID.
    std::vector<
        std::vector<std::reference_wrapper<MeshLib::PropertyVector<double>>>>
        residuum_vectors;

    virtual ~CommonAssemblyData() = default;

protected:
    bool areAllElementsActive() const
    {
        return sorted_active_element_ids_ == nullptr;
    }

    // shared_ptr allows to return either this data or data computed on the fly
    // from a member function (cf. activeElementIDsSorted()).
    std::shared_ptr<std::vector<std::size_t> const> sorted_active_element_ids_;
};

//! Data necessary for global equation system assembly on the bulk mesh.
struct BulkMeshAssemblyData : CommonAssemblyData
{
    using CommonAssemblyData::CommonAssemblyData;

    //! Assembly should proceed on all mesh elements.
    void setAllElementsActive();

    //! Assembly should proceed on the passed element IDs only.
    void setElementSelectionActive(
        std::vector<std::size_t> const& sorted_active_element_ids_whole_mesh);
};

//! Data necessary for global equation system assembly on submeshes of the bulk
//! mesh.
struct SubmeshAssemblyData : CommonAssemblyData
{
    //! Creates a new instance.
    //!
    //! \param submesh The submesh to which the instance belongs.
    //! \param residuum_vectors Residuum vectors for each process ID.
    explicit SubmeshAssemblyData(
        MeshLib::Mesh const& submesh,
        std::vector<std::vector<
            std::reference_wrapper<MeshLib::PropertyVector<double>>>>&&
            residuum_vectors);

    //! Assembly should proceed on all elements of this submesh.
    void setAllElementsActive();

    //! Assembly should proceed on the passed element IDs only (intersected with
    //! the elements of this submesh).
    void setElementSelectionActive(
        std::vector<std::size_t> const& sorted_active_element_ids_whole_mesh);

    MeshLib::PropertyVector<std::size_t> const& bulk_node_ids;

private:
    MeshLib::PropertyVector<std::size_t> const& bulk_element_ids;
};
}  // namespace ProcessLib::Assembly
