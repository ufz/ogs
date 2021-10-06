/*
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <algorithm>
#include <cstdlib>
#include <vector>

#include "GeoLib/GEOObjects.h"
#include "GeoLib/Polygon.h"
#include "MeshGeoToolsLib/MeshEditing/MarkNodesOutsideOfPolygon.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/IO/writeMeshToFile.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"

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

    auto const outside = markNodesOutSideOfPolygon(mesh.getNodes(), polygon);

    auto is_node_outside = [&outside](auto const* node_ptr)
    { return outside[node_ptr->getID()]; };

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

    using func =
        std::function<bool(MeshLib::Node* const*, MeshLib::Node* const*,
                           decltype(is_node_outside))>;
    auto is_element_outside =
        any_of ? (func)(std::all_of<MeshLib::Node* const*,
                                    decltype(is_node_outside)>)
               : (func)(std::any_of<MeshLib::Node* const*,
                                    decltype(is_node_outside)>);

    for (std::size_t j(0); j < mesh.getElements().size(); ++j)
    {
        MeshLib::Element const* const elem(mesh.getElements()[j]);
        if (is_element_outside(elem->getNodes(),
                               elem->getNodes() + elem->getNumberOfNodes(),
                               is_node_outside))
        {
            continue;
        }
        if (has_element_required_material_id(elem->getID()))
        {
            (*pv)[j] = new_property_value;
        }
    }
}

}  // namespace MeshGeoToolsLib
