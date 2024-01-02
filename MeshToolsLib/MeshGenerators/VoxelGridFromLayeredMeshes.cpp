/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "VoxelGridFromLayeredMeshes.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/MeshSearch/MeshElementGrid.h"
#include "MeshToolsLib/MeshEditing/ProjectPointOnMesh.h"
#include "MeshToolsLib/MeshEditing/RemoveMeshComponents.h"
#include "MeshToolsLib/MeshGenerators/MeshGenerator.h"

static std::string mat_name = "MaterialIDs";

// returns the AABB of all mesh nodes of layers read so far
void adjustExtent(std::pair<MathLib::Point3d, MathLib::Point3d>& extent,
                  MeshLib::Mesh const& mesh)
{
    auto const& nodes = mesh.getNodes();
    GeoLib::AABB aabb(nodes.cbegin(), nodes.cend());
    for (std::size_t i = 0; i < 3; ++i)
    {
        extent.first[i] = std::min(extent.first[i], aabb.getMinPoint()[i]);
        extent.second[i] = std::max(extent.second[i], aabb.getMaxPoint()[i]);
    }
}

// creates a voxel grid of the AABB of all layers
std::unique_ptr<MeshLib::Mesh> generateInitialMesh(
    std::pair<MathLib::Point3d, MathLib::Point3d>& extent,
    std::array<double, 3> const& res)
{
    INFO("Creating initial mesh...");
    std::array<double, 3> mesh_range{{extent.second[0] - extent.first[0],
                                      extent.second[1] - extent.first[1],
                                      extent.second[2] - extent.first[2]}};
    std::array<std::size_t, 3> const n_cells{
        {static_cast<std::size_t>(std::ceil(mesh_range[0] / res[0])),
         static_cast<std::size_t>(std::ceil(mesh_range[1] / res[1])),
         static_cast<std::size_t>(std::ceil(mesh_range[2] / res[2]))}};
    for (std::size_t i = 0; i < 3; ++i)
    {
        double const ext_range = n_cells[i] * res[i];
        double const offset = (ext_range - mesh_range[i]) / 2.0;
        mesh_range[i] = ext_range;
        extent.first[i] -= offset;
        extent.second[i] += offset;
    }
    std::unique_ptr<MeshLib::Mesh> mesh(
        MeshToolsLib::MeshGenerator::generateRegularHexMesh(
            mesh_range[0], mesh_range[1], mesh_range[2], n_cells[0], n_cells[1],
            n_cells[2], extent.first));
    auto mat_id = mesh->getProperties().createNewPropertyVector<int>(
        mat_name, MeshLib::MeshItemType::Cell);
    if (!mat_id)
    {
        return nullptr;
    }
    mat_id->insert(mat_id->end(), mesh->getNumberOfElements(), -1);
    return mesh;
}

// returns the element the given node is projected on (or nullptr otherwise)
MeshLib::Element const* getProjectedElement(
    MeshLib::MeshElementGrid const& grid,
    MathLib::Point3d const& node,
    double const max_edge)
{
    constexpr double max_val = std::numeric_limits<double>::max();
    MathLib::Point3d const min_vol{
        {node[0] - max_edge, node[1] - max_edge, -max_val}};
    MathLib::Point3d const max_vol{
        {node[0] + max_edge, node[1] + max_edge, max_val}};
    auto const& intersection_candidates =
        grid.getElementsInVolume(min_vol, max_vol);
    return MeshToolsLib::ProjectPointOnMesh::getProjectedElement(
        intersection_candidates, node);
}

// casts vote if the given nodes belongs to lower layer, upper layer or no layer
// at all
void voteMatId(MathLib::Point3d const& node,
               MeshLib::MeshElementGrid const& grid,
               double const max_edge,
               std::size_t& nullptr_cnt,
               std::size_t& upper_layer_cnt,
               std::size_t& lower_layer_cnt)
{
    auto const& proj_elem = getProjectedElement(grid, node, max_edge);
    if (proj_elem == nullptr)
    {
        nullptr_cnt++;
        return;
    }
    if (node[2] >
        MeshToolsLib::ProjectPointOnMesh::getElevation(*proj_elem, node))
    {
        upper_layer_cnt++;
        return;
    }
    lower_layer_cnt++;
}

// sets material IDs for all elements depending on the layers they are located
// between
void setMaterialIDs(MeshLib::Mesh& mesh,
                    std::vector<MeshLib::Mesh const*> const& layers,
                    bool const dilate)
{
    INFO("Setting material properties...");
    std::size_t const n_layers = layers.size();
    auto const& elems = mesh.getElements();
    std::size_t const n_elems = mesh.getNumberOfElements();
    auto mat_ids = mesh.getProperties().getPropertyVector<int>(mat_name);
    std::vector<bool> is_set(n_elems, false);
    for (int i = n_layers - 1; i >= 0; --i)
    {
        INFO("-> Layer {:d}", n_layers - i - 1);
        MeshLib::MeshElementGrid const grid(*layers[i]);
        auto const edgeLengths = minMaxEdgeLength(layers[i]->getElements());
        double const max_edge = edgeLengths.second;
        for (std::size_t j = 0; j < n_elems; ++j)
        {
            if (is_set[j])
            {
                continue;
            }

            std::size_t nullptr_cnt(0);
            std::size_t upper_layer_cnt(0);
            std::size_t lower_layer_cnt(0);

            auto const& node = MeshLib::getCenterOfGravity(*elems[j]);
            voteMatId(node, grid, max_edge, nullptr_cnt, upper_layer_cnt,
                      lower_layer_cnt);
            if (nullptr_cnt)
            {
                // if no element was found at centre point, vote via corners
                for (std::size_t k = 0; k < 8; ++k)
                {
                    MeshLib::Node const& n = *elems[j]->getNode(k);
                    voteMatId(n, grid, max_edge, nullptr_cnt, upper_layer_cnt,
                              lower_layer_cnt);
                }

                // If the "dilate"-param is set, a mat ID will be assigned if at
                // least one node was voting for a specific layer. Without the
                // "dilate"-param, an absolute majority is needed. In case of a
                // tie, the lower layer will be favoured.
                if ((upper_layer_cnt == 0 && lower_layer_cnt == 0) ||
                    (!dilate && nullptr_cnt >= upper_layer_cnt &&
                     nullptr_cnt >= lower_layer_cnt))
                {
                    continue;
                }
                if (upper_layer_cnt > lower_layer_cnt)
                {
                    (*mat_ids)[j] = n_layers - i - 1;
                }
                else
                {
                    is_set[j] = true;
                }
                continue;
            }
            if (upper_layer_cnt)
            {
                (*mat_ids)[j] = n_layers - i - 1;
            }
            else
            {
                is_set[j] = true;
            }
        }
    }
    // set all elements above uppermost layer back to -1 so they are
    // subsequently cut
    std::replace(mat_ids->begin(), mat_ids->end(),
                 static_cast<int>(n_layers - 1), -1);
}

// Removes all elements from mesh that have not been marked as being located
// between two layers. If all elements remain unmarked, a nullptr is returned.
std::vector<std::size_t> markSpecificElements(MeshLib::Mesh const& mesh,
                                              int const mat_id)
{
    std::vector<std::size_t> marked_elems;
    auto const mat_ids = *MeshLib::materialIDs(mesh);
    std::size_t const n_elems = mat_ids.size();
    for (std::size_t i = 0; i < n_elems; ++i)
    {
        if (mat_ids[i] == mat_id)
        {
            marked_elems.push_back(i);
        }
    }
    return marked_elems;
}

// Creates a VoxelGrid after extending the AABB for each layer.
std::unique_ptr<MeshLib::Mesh> MeshToolsLib::MeshGenerators::
    VoxelFromLayeredMeshes::createVoxelFromLayeredMesh(
        std::pair<MathLib::Point3d, MathLib::Point3d>& extent,
        std::vector<MeshLib::Mesh const*> const& layers,
        std::array<double, 3> const cellsize,
        bool const dilate)
{
    for (auto const& layer : layers)
    {
        adjustExtent(extent, *layer);
    }

    std::unique_ptr<MeshLib::Mesh> mesh(generateInitialMesh(extent, cellsize));
    if (mesh == nullptr)
    {
        return nullptr;
    }
    setMaterialIDs(*mesh, layers, dilate);
    auto const marked_elements = markSpecificElements(*mesh, -1);
    if (marked_elements.size() == mesh->getNumberOfElements())
    {
        return nullptr;
    }
    std::unique_ptr<MeshLib::Mesh> new_mesh(
        MeshToolsLib::removeElements(*mesh, marked_elements, "mesh"));
    return new_mesh;
}
