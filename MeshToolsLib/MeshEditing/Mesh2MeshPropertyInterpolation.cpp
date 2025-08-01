/**
 * \file
 * \author Thomas Fischer
 * \date   Oct 12, 2012
 * \brief  Implementation of the Mesh2MeshPropertyInterpolation class.
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Mesh2MeshPropertyInterpolation.h"

#include <fstream>
#include <vector>

#include "BaseLib/Logging.h"
#include "GeoLib/AABB.h"
#include "GeoLib/Grid.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshEnums.h"
#include "MeshLib/Node.h"

namespace MeshToolsLib
{
Mesh2MeshPropertyInterpolation::Mesh2MeshPropertyInterpolation(
    MeshLib::Mesh const& src_mesh, std::string const& property_name)
    : _src_mesh(src_mesh), _property_name(property_name)
{
}

bool Mesh2MeshPropertyInterpolation::setPropertiesForMesh(
    MeshLib::Mesh& dest_mesh) const
{
    if (_src_mesh.getDimension() != dest_mesh.getDimension())
    {
        ERR("MeshLib::Mesh2MeshPropertyInterpolation::setPropertiesForMesh() "
            "dimension of source (dim = {:d}) and destination (dim = {:d}) "
            "mesh does not match.",
            _src_mesh.getDimension(), dest_mesh.getDimension());
        return false;
    }

    if (_src_mesh.getDimension() != 2)
    {
        WARN(
            "MeshLib::Mesh2MeshPropertyInterpolation::setPropertiesForMesh() "
            "implemented only for 2D case at the moment.");
        return false;
    }

    MeshLib::PropertyVector<double>* dest_properties;
    if (dest_mesh.getProperties().existsPropertyVector<double>(_property_name))
    {
        dest_properties =
            dest_mesh.getProperties().getPropertyVector<double>(_property_name);
        assert(dest_properties->size() != dest_mesh.getNumberOfElements());
    }
    else
    {
        INFO("Create new PropertyVector '{:s}' of type double.",
             _property_name);
        dest_properties =
            dest_mesh.getProperties().createNewPropertyVector<double>(
                _property_name, MeshLib::MeshItemType::Cell,
                dest_mesh.getNumberOfElements(), 1);
        if (!dest_properties)
        {
            WARN(
                "Could not get or create a PropertyVector of type double using "
                "the given name '{:s}'.",
                _property_name);
            return false;
        }
    }

    interpolatePropertiesForMesh(dest_mesh, *dest_properties);

    return true;
}

void Mesh2MeshPropertyInterpolation::interpolatePropertiesForMesh(
    MeshLib::Mesh const& dest_mesh,
    MeshLib::PropertyVector<double>& dest_properties) const
{
    std::vector<double> interpolated_src_node_properties(
        _src_mesh.getNumberOfNodes());
    interpolateElementPropertiesToNodeProperties(
        interpolated_src_node_properties);

    // idea: looping over the destination elements and calculate properties
    // from interpolated_src_node_properties to accelerate the (source) point
    // search construct a grid
    std::vector<MeshLib::Node*> const& src_nodes(_src_mesh.getNodes());
    GeoLib::Grid<MeshLib::Node> src_grid(src_nodes.begin(), src_nodes.end(),
                                         64);

    auto const& dest_elements(dest_mesh.getElements());
    const std::size_t n_dest_elements(dest_elements.size());
    for (std::size_t k(0); k < n_dest_elements; k++)
    {
        MeshLib::Element& dest_element(*dest_elements[k]);
        if (dest_element.getGeomType() == MeshLib::MeshElemType::LINE)
        {
            continue;
        }

        // compute axis aligned bounding box around the current element
        const GeoLib::AABB elem_aabb(
            dest_element.getNodes(),
            dest_element.getNodes() + dest_element.getNumberOfBaseNodes());

        // request "interesting" nodes from grid
        std::vector<std::vector<MeshLib::Node*> const*> const nodes =
            src_grid.getPntVecsOfGridCellsIntersectingCuboid(
                elem_aabb.getMinPoint(), elem_aabb.getMaxPoint());

        std::size_t cnt(0);
        double average_value(0.0);

        for (auto const* nodes_vec : nodes)
        {
            for (auto const* node : *nodes_vec)
            {
                if (elem_aabb.containsPointXY(*node) &&
                    MeshLib::isPointInElementXY(*node, dest_element))
                {
                    average_value +=
                        interpolated_src_node_properties[node->getID()];
                    cnt++;
                }
            }
        }

        if (cnt == 0)
        {
            OGS_FATAL(
                "Mesh2MeshInterpolation: Could not find values in source mesh "
                "for the element {:d}.",
                k);
        }
        dest_properties[k] = average_value / cnt;
    }
}

void Mesh2MeshPropertyInterpolation::
    interpolateElementPropertiesToNodeProperties(
        std::vector<double>& interpolated_properties) const
{
    // fetch the source of property values
    if (!_src_mesh.getProperties().existsPropertyVector<double>(_property_name))
    {
        WARN("Did not find PropertyVector<double> '{:s}'.", _property_name);
        return;
    }
    auto const* elem_props =
        _src_mesh.getProperties().getPropertyVector<double>(_property_name);

    std::vector<MeshLib::Node*> const& src_nodes(_src_mesh.getNodes());
    const std::size_t n_src_nodes(src_nodes.size());
    for (std::size_t k(0); k < n_src_nodes; k++)
    {
        const std::size_t n_con_elems(
            _src_mesh.getElementsConnectedToNode(*src_nodes[k]).size());
        interpolated_properties[k] = (*elem_props)
            [_src_mesh.getElementsConnectedToNode(*src_nodes[k])[0]->getID()];
        for (std::size_t j(1); j < n_con_elems; j++)
        {
            interpolated_properties[k] +=
                (*elem_props)[_src_mesh
                                  .getElementsConnectedToNode(*src_nodes[k])[j]
                                  ->getID()];
        }
        interpolated_properties[k] /= n_con_elems;
    }
}

}  // namespace MeshToolsLib
