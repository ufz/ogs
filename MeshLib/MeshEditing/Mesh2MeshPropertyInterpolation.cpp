/**
 * \file
 * \author Thomas Fischer
 * \date   Oct 12, 2012
 * \brief  Implementation of the Mesh2MeshPropertyInterpolation class.
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <vector>
#include <fstream>
#include <boost/optional.hpp>

#include "Mesh2MeshPropertyInterpolation.h"

#include <logog/include/logog.hpp>

#include "GeoLib/AABB.h"
#include "GeoLib/Grid.h"

#include "MeshLib/MeshEnums.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/Elements/Element.h"

namespace MeshLib {

Mesh2MeshPropertyInterpolation::Mesh2MeshPropertyInterpolation(
    Mesh const*const src_mesh, std::string const& property_name)
    : _src_mesh(src_mesh), _property_name(property_name)
{}

bool Mesh2MeshPropertyInterpolation::setPropertiesForMesh(Mesh *dest_mesh) const
{
    if (_src_mesh->getDimension() != dest_mesh->getDimension()) {
        ERR ("MeshLib::Mesh2MeshPropertyInterpolation::setPropertiesForMesh() dimension of source (dim = %d) and destination (dim = %d) mesh does not match.", _src_mesh->getDimension(), dest_mesh->getDimension());
        return false;
    }

    if (_src_mesh->getDimension() != 2) {
        WARN ("MeshLib::Mesh2MeshPropertyInterpolation::setPropertiesForMesh() implemented only for 2D case at the moment.");
        return false;
    }

    GeoLib::AABB src_aabb(_src_mesh->getNodes().begin(), _src_mesh->getNodes().end());
    GeoLib::AABB dest_aabb(dest_mesh->getNodes().begin(), dest_mesh->getNodes().end());
    if (!src_aabb.containsAABB(dest_aabb)) {
        ERR("MeshLib::Mesh2MeshPropertyInterpolation::setPropertiesForMesh() source mesh to small.");
        ERR("src_aabb: %f, %f, %f | %f, %f, %f", src_aabb.getMinPoint()[0], src_aabb.getMinPoint()[1], src_aabb.getMinPoint()[2], src_aabb.getMaxPoint()[0], src_aabb.getMaxPoint()[1], src_aabb.getMaxPoint()[2]);
        ERR("dest_aabb: %f, %f, %f | %f, %f, %f", dest_aabb.getMinPoint()[0], dest_aabb.getMinPoint()[1], dest_aabb.getMinPoint()[2], dest_aabb.getMaxPoint()[0], dest_aabb.getMaxPoint()[1], dest_aabb.getMaxPoint()[2]);
        return false;
    }

    interpolatePropertiesForMesh(dest_mesh);

    return true;
}

void Mesh2MeshPropertyInterpolation::interpolatePropertiesForMesh(
    Mesh *dest_mesh) const
{
    // check the existence of PropertyVector in the source mesh
    boost::optional<MeshLib::PropertyVector<double> const&> opt_src_pv(
        _src_mesh->getProperties().getPropertyVector<double>(_property_name));
    if (!opt_src_pv) {
        WARN("Did not find PropertyVector<double> \"%s\" in source mesh.",
            _property_name.c_str());
        return;
    }

    // carry over property information from source elements to source nodes
    std::vector<double> interpolated_src_node_properties(_src_mesh->getNumberOfNodes());
    interpolateElementPropertiesToNodeProperties(interpolated_src_node_properties);

    boost::optional<MeshLib::PropertyVector<double> &> opt_pv(
        dest_mesh->getProperties().getPropertyVector<double>(_property_name));
    if (!opt_pv) {
        INFO("Create new PropertyVector \"%s\" of type double.",
             _property_name.c_str());
        opt_pv = dest_mesh->getProperties().createNewPropertyVector<double>(
            _property_name, MeshItemType::Cell, 1);
        if (!opt_pv) {
            WARN("Could not get or create a PropertyVector of type double"
                " using the given name \"%s\".", _property_name.c_str());
            return;
        }
    }
    MeshLib::PropertyVector<double> & dest_properties(opt_pv.get());
    if (dest_properties.size() != dest_mesh->getNumberOfElements())
        dest_properties.resize(dest_mesh->getNumberOfElements());

    // looping over the destination elements and calculate properties
    // from interpolated_src_node_properties
    std::vector<MeshLib::Node*> const& src_nodes(_src_mesh->getNodes());

    GeoLib::Grid<MeshLib::Node> src_grid(src_nodes.begin(), src_nodes.end(), 64);

    std::vector<MeshLib::Element*> const& dest_elements(dest_mesh->getElements());
    const std::size_t n_dest_elements(dest_elements.size());
    for (std::size_t k(0); k<n_dest_elements; k++)
    {
        if (dest_elements[k]->getGeomType() == MeshElemType::LINE)
            continue;

        // compute axis aligned bounding box around the current element
        const GeoLib::AABB elem_aabb(
            dest_elements[k]->getNodes(),
            dest_elements[k]->getNodes() +
                dest_elements[k]->getNumberOfBaseNodes());

        // request "interesting" nodes from grid
        std::vector<std::vector<MeshLib::Node*> const*> nodes;
        src_grid.getPntVecsOfGridCellsIntersectingCuboid(
            elem_aabb.getMinPoint(), elem_aabb.getMaxPoint(), nodes);

        std::size_t cnt(0);
        dest_properties[k] = 0.0;

        for (auto i_th_vec : nodes) {
            const std::size_t n_nodes_in_vec(i_th_vec->size());
            for (std::size_t j(0); j<n_nodes_in_vec; j++) {
                MeshLib::Node const*const j_th_node((*i_th_vec)[j]);
                if (elem_aabb.containsPoint(*j_th_node)) {
                    if (dest_elements[k]->isPntInElement(*j_th_node)) {
                        dest_properties[k] += interpolated_src_node_properties[(*i_th_vec)[j]->getID()];
                        cnt++;
                    }
                }
            }
        }

        dest_properties[k] /= cnt;

        if (cnt == 0) {
            std::string base_name("DebugData/Element-");
            base_name += std::to_string(k);

            std::string aabb_fname(base_name + "-aabb.gli");
            std::ofstream out_aabb(aabb_fname.c_str());
            out_aabb << "#POINTS" << "\n";
            out_aabb << "0 " << elem_aabb.getMinPoint() << "\n";
            out_aabb << "1 " << elem_aabb.getMaxPoint() << "\n";
            out_aabb << "#STOP" << "\n";
            out_aabb.close();


            std::string source_fname(base_name + "-SourceNodes.gli");
            std::ofstream out_src(source_fname.c_str());
            out_src << "#POINTS" << "\n";
            std::size_t nodes_cnt(0);
            for (auto i_th_vec : nodes) {
                const std::size_t n_nodes_in_vec(i_th_vec->size());
                for (std::size_t j(0); j<n_nodes_in_vec; j++) {
                    MeshLib::Node const*const j_th_node((*i_th_vec)[j]);
                    out_src << nodes_cnt << " " << *j_th_node << "\n";
                    nodes_cnt++;
                }
            }
            out_src << "#STOP" << "\n";
            out_src.close();
            ERR("no source nodes in dest element %d", k);
        }
    }
}

void Mesh2MeshPropertyInterpolation::interpolateElementPropertiesToNodeProperties(
    std::vector<double> &interpolated_properties) const
{
    // fetch the source of property values
    boost::optional<MeshLib::PropertyVector<double> const&> opt_src_props(
        _src_mesh->getProperties().getPropertyVector<double>(_property_name));
    if (!opt_src_props)
    {
        WARN("Did not find PropertyVector<double> \"%s\".",
             _property_name.c_str());
        return;
    }

    MeshLib::PropertyVector<double> const& elem_props(opt_src_props.get());
    std::vector<MeshLib::Node*> const& src_nodes(_src_mesh->getNodes());
    const std::size_t n_src_nodes(src_nodes.size());
    for (std::size_t k(0); k < n_src_nodes; k++)
    {
        const std::size_t n_con_elems(src_nodes[k]->getNumberOfElements());
        interpolated_properties[k] = elem_props[(src_nodes[k]->getElement(0))->getID()];
        for (std::size_t j(1); j < n_con_elems; j++)
        {
            interpolated_properties[k] += elem_props[(src_nodes[k]->getElement(j))->getID()];
        }
        interpolated_properties[k] /= n_con_elems;
    }
}

} // end namespace MeshLib
