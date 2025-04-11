/*
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <algorithm>
#include <cstdlib>
#include <range/v3/algorithm/all_of.hpp>
#include <range/v3/algorithm/any_of.hpp>
#include <range/v3/view/transform.hpp>
#include <vector>

#include "GeoLib/GEOObjects.h"
#include "GeoLib/Polygon.h"
#include "MeshGeoToolsLib/MeshEditing/MarkNodesOutsideOfPolygon.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/IO/writeMeshToFile.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/Utils/getOrCreateMeshProperty.h"

namespace MeshGeoToolsLib
{
template <typename PT>
void resetMeshElementProperty(MeshLib::Mesh& mesh,
                              GeoLib::Polygon const& polygon,
                              std::string const& property_name,
                              PT new_property_value,
                              int restrict_to_material_id,
                              bool const any_of)
{
    auto* const pv = MeshLib::getOrCreateMeshProperty<PT>(
        mesh, property_name, MeshLib::MeshItemType::Cell, 1);

    if (pv->getMeshItemType() != MeshLib::MeshItemType::Cell)
    {
        ERR("Values of the PropertyVector are not assigned to cells.");
        return;
    }
    auto const* material_ids =
        mesh.getProperties().getPropertyVector<int>("MaterialIDs");

    if (restrict_to_material_id != -1 && !material_ids)
    {
        OGS_FATAL(
            "Restriction of resetting a property in a polygonal region "
            "requires that a MaterialIDs data array is available in the "
            "mesh.");
    }

    auto has_element_required_material_id = [&](int const element_id)
    {
        return restrict_to_material_id == -1 ||
               (*material_ids)[element_id] == restrict_to_material_id;
    };

    auto const outside = markNodesOutSideOfPolygon(mesh.getNodes(), polygon);

    auto is_node_outside = [&outside](std::size_t const node_id)
    { return outside[node_id]; };

    auto is_element_outside = [&](MeshLib::Element const& e)
    {
        auto ids = e.nodes() | MeshLib::views::ids;
        return any_of ? ranges::all_of(ids, is_node_outside)
                      : ranges::any_of(ids, is_node_outside);
    };

    auto is_valid_element = [&](MeshLib::Element const& e)
    {
        return !is_element_outside(e) &&
               has_element_required_material_id(e.getID());
    };

    auto compute_value = [&](MeshLib::Element const* const e) -> PT
    { return is_valid_element(*e) ? new_property_value : (*pv)[e->getID()]; };

    pv->assign(mesh.getElements() | ranges::views::transform(compute_value));
}

}  // namespace MeshGeoToolsLib
