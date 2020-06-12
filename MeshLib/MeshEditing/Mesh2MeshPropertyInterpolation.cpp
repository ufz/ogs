/**
 * \file
 * \author Thomas Fischer
 * \date   Oct 12, 2012
 * \brief  Implementation of the Mesh2MeshPropertyInterpolation class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <vector>
#include <fstream>

#include "Mesh2MeshPropertyInterpolation.h"

#include "BaseLib/Logging.h"

#include "GeoLib/AABB.h"
#include "GeoLib/Grid.h"

#include "MeshLib/MeshEnums.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/Elements/Element.h"

namespace MeshLib {

Mesh2MeshPropertyInterpolation::Mesh2MeshPropertyInterpolation(
    Mesh const& src_mesh, std::string const& property_name)
    : src_mesh_(src_mesh), property_name_(property_name)
{}

bool Mesh2MeshPropertyInterpolation::setPropertiesForMesh(Mesh& dest_mesh) const
{
    if (src_mesh_.getDimension() != dest_mesh.getDimension()) {
        ERR("MeshLib::Mesh2MeshPropertyInterpolation::setPropertiesForMesh() "
            "dimension of source (dim = {:d}) and destination (dim = {:d}) "
            "mesh does not match.",
            src_mesh_.getDimension(), dest_mesh.getDimension());
        return false;
    }

    if (src_mesh_.getDimension() != 2) {
        WARN(
            "MeshLib::Mesh2MeshPropertyInterpolation::setPropertiesForMesh() "
            "implemented only for 2D case at the moment.");
        return false;
    }

    MeshLib::PropertyVector<double>* dest_properties;
    if (dest_mesh.getProperties().existsPropertyVector<double>(property_name_))
    {
        dest_properties =
            dest_mesh.getProperties().getPropertyVector<double>(property_name_);
    }
    else
    {
        INFO("Create new PropertyVector '{:s}' of type double.",
             property_name_);
        dest_properties =
            dest_mesh.getProperties().createNewPropertyVector<double>(
                property_name_, MeshItemType::Cell, 1);
        if (!dest_properties)
        {
            WARN(
                "Could not get or create a PropertyVector of type double"
                " using the given name '{:s}'.",
                property_name_);
            return false;
        }
    }
    if (dest_properties->size() != dest_mesh.getNumberOfElements())
    {
        dest_properties->resize(dest_mesh.getNumberOfElements());
    }

    interpolatePropertiesForMesh(dest_mesh, *dest_properties);

    return true;
}

void Mesh2MeshPropertyInterpolation::interpolatePropertiesForMesh(
    Mesh const& dest_mesh,
    MeshLib::PropertyVector<double>& dest_properties) const
{
    std::vector<double> interpolated_src_node_properties(
        src_mesh_.getNumberOfNodes());
    interpolateElementPropertiesToNodeProperties(
        interpolated_src_node_properties);

    // idea: looping over the destination elements and calculate properties
    // from interpolated_src_node_properties to accelerate the (source) point
    // search construct a grid
    std::vector<MeshLib::Node*> const& src_nodes(src_mesh_.getNodes());
    GeoLib::Grid<MeshLib::Node> src_grid(src_nodes.begin(), src_nodes.end(),
                                         64);

    auto const& dest_elements(dest_mesh.getElements());
    const std::size_t n_dest_elements(dest_elements.size());
    for (std::size_t k(0); k < n_dest_elements; k++)
    {
        MeshLib::Element& dest_element(*dest_elements[k]);
        if (dest_element.getGeomType() == MeshElemType::LINE)
        {
            continue;
        }

        // compute axis aligned bounding box around the current element
        const GeoLib::AABB elem_aabb(
            dest_element.getNodes(),
            dest_element.getNodes() + dest_element.getNumberOfBaseNodes());

        // request "interesting" nodes from grid
        std::vector<std::vector<MeshLib::Node*> const*> nodes;
        src_grid.getPntVecsOfGridCellsIntersectingCuboid(
            elem_aabb.getMinPoint(), elem_aabb.getMaxPoint(), nodes);

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

void Mesh2MeshPropertyInterpolation::interpolateElementPropertiesToNodeProperties(
    std::vector<double> &interpolated_properties) const
{
    // fetch the source of property values
    if (!src_mesh_.getProperties().existsPropertyVector<double>(property_name_))
    {
        WARN("Did not find PropertyVector<double> '{:s}'.", property_name_);
        return;
    }
    auto const* elem_props =
        src_mesh_.getProperties().getPropertyVector<double>(property_name_);

    std::vector<MeshLib::Node*> const& src_nodes(src_mesh_.getNodes());
    const std::size_t n_src_nodes(src_nodes.size());
    for (std::size_t k(0); k < n_src_nodes; k++)
    {
        const std::size_t n_con_elems(src_nodes[k]->getNumberOfElements());
        interpolated_properties[k] =
            (*elem_props)[(src_nodes[k]->getElement(0))->getID()];
        for (std::size_t j(1); j < n_con_elems; j++)
        {
            interpolated_properties[k] +=
                (*elem_props)[(src_nodes[k]->getElement(j))->getID()];
        }
        interpolated_properties[k] /= n_con_elems;
    }
}

} // end namespace MeshLib
