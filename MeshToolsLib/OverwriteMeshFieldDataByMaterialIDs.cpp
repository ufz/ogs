// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "OverwriteMeshFieldDataByMaterialIDs.h"

#include <cstddef>
#include <vector>

#include "BaseLib/Error.h"
#include "BaseLib/Logging.h"
#include "IntegrationPointDataTools.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Utils/IntegrationPointWriter.h"
#include "ParameterLib/SpatialPosition.h"

namespace MeshToolsLib
{
namespace
{
//! A contiguous run of "slots" (the integration points of one element, a single
//! cell, or a single node) that all receive the same value when a parameter is
//! used. \c slot multiplied by the number of components gives the flat offset
//! into the property vector.
struct WriteGroup
{
    //! Position at which a parameter is evaluated for this group.
    ParameterLib::SpatialPosition pos;
    std::size_t slot_begin;
    std::size_t slot_end;  //!< exclusive
};

//! Collects the write groups for \c initial_condition according to its mesh
//! item type. The three supported item types (integration point, cell, node)
//! differ only in which slots are written and where a parameter is evaluated;
//! everything else is shared by overwriteMeshFieldDataByMaterialIDs().
std::vector<WriteGroup> collectWriteGroups(
    InitialConditionDataSet const& initial_condition,
    MeshLib::PropertyVector<double> const& pv,
    MeshLib::Properties const& properties)
{
    auto const& mesh = *initial_condition.mesh;
    auto const& element_ids =
        initial_condition.element_ids_for_selected_materials;

    std::vector<WriteGroup> groups;

    switch (initial_condition.mesh_item_type)
    {
        case MeshLib::MeshItemType::IntegrationPoint:
        {
            auto const element_ip_data_offsets =
                getIntegrationPointDataOffsetsOfMeshElements(mesh.getElements(),
                                                             pv, properties);
            // getIntegrationPointDataOffsetsOfMeshElements() returns an empty
            // vector for a property whose name does not contain "_ip" (it is
            // then not treated as integration-point data). Indexing it by
            // element id below would be an out-of-bounds access.
            if (element_ip_data_offsets.empty())
            {
                OGS_FATAL(
                    "overwrite_mesh_data: integration-point property '{:s}' "
                    "carries no integration-point data offsets. Only fields "
                    "whose name contains '_ip' are treated as "
                    "integration-point data.",
                    initial_condition.variable_name);
            }
            auto const n = pv.getNumberOfGlobalComponents();
            groups.reserve(element_ids.size());
            for (auto const element_id : element_ids)
            {
                // A parameter is evaluated once per element at the centroid and
                // written to all of the element's integration points (a single
                // write group spanning them). This is exact for the intended
                // constant-per-material use case. Do not turn this into a
                // per-integration-point evaluation without also reproducing the
                // exact integration rule and point ordering used by the process
                // local assembler that wrote this field: computing the points
                // in a different order would associate values with the wrong
                // slots.
                ParameterLib::SpatialPosition pos;
                pos.setElementID(element_id);
                pos.setCoordinates(
                    MeshLib::getCenterOfGravity(*mesh.getElement(element_id)));
                // The offsets are flat (already multiplied by the number of
                // components); divide to obtain slot (integration point)
                // indices.
                groups.push_back({pos, element_ip_data_offsets[element_id] / n,
                                  element_ip_data_offsets[element_id + 1] / n});
            }
            break;
        }
        case MeshLib::MeshItemType::Cell:
        {
            groups.reserve(element_ids.size());
            for (auto const element_id : element_ids)
            {
                ParameterLib::SpatialPosition pos;
                pos.setElementID(element_id);
                pos.setCoordinates(
                    MeshLib::getCenterOfGravity(*mesh.getElement(element_id)));
                groups.push_back({pos, element_id, element_id + 1});
            }
            break;
        }
        case MeshLib::MeshItemType::Node:
        {
            for (auto const element_id : element_ids)
            {
                for (auto const* node : mesh.getElement(element_id)->nodes())
                {
                    std::size_t const node_id = node->getID();
                    // Evaluate the parameter at the node position, not the
                    // element centre.
                    groups.push_back(
                        {ParameterLib::SpatialPosition{node_id, {}, *node},
                         node_id, node_id + 1});
                }
            }
            break;
        }
        default:
            OGS_FATAL(
                "The mesh item type of the initial condition for property "
                "'{:s}' is not supported. Only integration point, cell, and "
                "node data are supported.",
                initial_condition.variable_name);
    }

    return groups;
}
}  // namespace

void overwriteMeshFieldDataByMaterialIDs(
    std::vector<InitialConditionDataSet>& data_initial_conditions)
{
    // First pass: snapshot the original field data for every non-parameter
    // (i.e. take_original) condition before any writes happen, so that a later
    // take_original can restore the original values even if an earlier set
    // overwrote overlapping slots.
    for (auto& initial_condition : data_initial_conditions)
    {
        auto const& property_name = initial_condition.variable_name;
        MeshLib::Properties& properties =
            initial_condition.mesh->getProperties();

        if (!properties.hasPropertyVector<double>(
                property_name, initial_condition.mesh_item_type))
        {
            OGS_FATAL(
                "overwrite_mesh_data: property '{:s}' of the requested mesh "
                "item type was not found on the mesh.",
                property_name);
        }

        // With a parameter the values are computed, so no snapshot is needed.
        if (initial_condition.parameter)
        {
            continue;
        }

        if (!initial_condition.property_copy)
        {
            auto* pv = properties.getPropertyVector<double>(property_name);
            initial_condition.property_copy.reset(
                static_cast<MeshLib::PropertyVector<double> const*>(
                    pv->clone({})));
            DBUG("created property copy for: {}", property_name);
        }
    }

    // Second pass: overwrite the selected slots for each condition.
    for (auto const& initial_condition : data_initial_conditions)
    {
        MeshLib::Properties& properties =
            initial_condition.mesh->getProperties();
        auto& pv = *properties.getPropertyVector<double>(
            initial_condition.variable_name);
        auto const n_components = pv.getNumberOfGlobalComponents();

        auto const groups =
            collectWriteGroups(initial_condition, pv, properties);

        bool const has_parameter = initial_condition.parameter != nullptr;
        // Either a parameter provides the values or a snapshot was taken in the
        // first pass; one of the two must hold. Without a snapshot copy_data
        // below is null and would be dereferenced.
        if (!has_parameter && !initial_condition.property_copy)
        {
            OGS_FATAL(
                "overwrite_mesh_data: neither a parameter nor an original-data "
                "snapshot is available for property '{:s}'. This is an "
                "internal "
                "inconsistency.",
                initial_condition.variable_name);
        }
        auto const* const copy_data =
            has_parameter ? nullptr
                          : static_cast<double const*>(
                                initial_condition.property_copy->data());

        for (auto const& group : groups)
        {
            std::vector<double> param_values;
            if (has_parameter)
            {
                param_values = (*initial_condition.parameter)(0.0, group.pos);
                // The parameter is looked up with num_components == 0 (the
                // component count is only known here), so validate it now:
                // reading n_components values from a shorter parameter result
                // would be an out-of-bounds access.
                if (static_cast<int>(param_values.size()) != n_components)
                {
                    OGS_FATAL(
                        "overwrite_mesh_data: the parameter used to set "
                        "property '{:s}' returns {:d} component(s), but the "
                        "property has {:d} component(s).",
                        initial_condition.variable_name, param_values.size(),
                        n_components);
                }
            }

            for (std::size_t slot = group.slot_begin; slot < group.slot_end;
                 ++slot)
            {
                double const* const source =
                    has_parameter ? param_values.data()
                                  : copy_data + slot * n_components;
                for (int comp = 0; comp < n_components; ++comp)
                {
                    pv[slot * n_components + comp] = source[comp];
                }
            }
        }
    }
}
}  // namespace MeshToolsLib
