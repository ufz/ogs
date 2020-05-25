/**
 * \file
 * \author Karsten Rink
 * \date   2014-04-11
 * \brief  Implementation of the LayeredVolume class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LayeredVolume.h"

#include "MathLib/Vector3.h"

#include "GeoLib/Raster.h"

#include "MeshLib/Elements/Tri.h"
#include "MeshLib/Elements/Quad.h"
#include "MeshLib/Properties.h"
#include "MeshLib/MeshEditing/DuplicateMeshComponents.h"
#include "MeshLib/MeshEditing/RemoveMeshComponents.h"
#include "MeshLib/MeshGenerators/MeshLayerMapper.h"
#include "MeshLib/MeshSearch/ElementSearch.h"


bool LayeredVolume::createRasterLayers(const MeshLib::Mesh &mesh,
                                       const std::vector<GeoLib::Raster const*> &rasters,
                                       double minimum_thickness,
                                       double noDataReplacementValue)
{
    if (mesh.getDimension() != 2)
    {
        return false;
    }

    elevation_epsilon_ = calcEpsilon(*rasters[0], *rasters.back());
    if (elevation_epsilon_ <= 0)
    {
        return false;
    }

    // remove line elements, only tri + quad remain
    MeshLib::ElementSearch ex(mesh);
    ex.searchByElementType(MeshLib::MeshElemType::LINE);
    std::unique_ptr<MeshLib::Mesh> top(
        removeElements(mesh, ex.getSearchedElementIDs(), "MeshLayer"));
    if (top == nullptr)
    {
        top = std::make_unique<MeshLib::Mesh>(mesh);
    }

    if (!MeshLib::MeshLayerMapper::layerMapping(
            *top, *rasters.back(), noDataReplacementValue))
    {
        return false;
    }

    std::unique_ptr<MeshLib::Mesh> bottom(new MeshLib::Mesh(*top));
    if (!MeshLib::MeshLayerMapper::layerMapping(*bottom, *rasters[0], 0))
    {
        return false;
    }

    this->minimum_thickness_ = minimum_thickness;
    nodes_ = MeshLib::copyNodeVector(bottom->getNodes());
    elements_ = MeshLib::copyElementVector(bottom->getElements(), nodes_);
    if (!materials_.empty())
    {
        ERR("The materials vector is not empty.");
        return false;
    }
    materials_.resize(elements_.size(), 0);

    // map each layer and attach to subsurface mesh
    const std::size_t nRasters (rasters.size());
    for (std::size_t i = 1; i < nRasters; ++i)
    {
        this->addLayerToMesh(*top, i, *rasters[i]);
    }

    // close boundaries between layers
    this->addLayerBoundaries(*top, nRasters);
    this->removeCongruentElements(nRasters, top->getNumberOfElements());

    return true;
}

void LayeredVolume::addLayerToMesh(const MeshLib::Mesh &dem_mesh, unsigned layer_id, GeoLib::Raster const& raster)
{
    const std::size_t nNodes (dem_mesh.getNumberOfNodes());
    const std::vector<MeshLib::Node*> &nodes (dem_mesh.getNodes());
    const std::size_t node_id_offset (nodes_.size());
    const std::size_t last_layer_node_offset (node_id_offset-nNodes);

    for (std::size_t i = 0; i < nNodes; ++i)
    {
        nodes_.push_back(getNewLayerNode(*nodes[i],
                                         *nodes_[last_layer_node_offset + i],
                                         raster,
                                         nodes_.size()));
    }

    const std::vector<MeshLib::Element*> &layer_elements (dem_mesh.getElements());
    for (MeshLib::Element* elem : layer_elements)
    {
        if (elem->getGeomType() == MeshLib::MeshElemType::TRIANGLE)
        {
            std::array<MeshLib::Node*,3> tri_nodes = {{ nodes_[node_id_offset+elem->getNodeIndex(0)],
                                                        nodes_[node_id_offset+elem->getNodeIndex(1)],
                                                        nodes_[node_id_offset+elem->getNodeIndex(2)] }};
            elements_.push_back(new MeshLib::Tri(tri_nodes));
            materials_.push_back(layer_id);
        }
        else if (elem->getGeomType() == MeshLib::MeshElemType::QUAD)
        {
            std::array<MeshLib::Node*,4> quad_nodes = {{ nodes_[node_id_offset+elem->getNodeIndex(0)],
                                                         nodes_[node_id_offset+elem->getNodeIndex(1)],
                                                         nodes_[node_id_offset+elem->getNodeIndex(2)],
                                                         nodes_[node_id_offset+elem->getNodeIndex(3)] }};
            elements_.push_back(new MeshLib::Quad(quad_nodes));
            materials_.push_back(layer_id);
        }
    }
}

void LayeredVolume::addLayerBoundaries(const MeshLib::Mesh &layer, std::size_t nLayers)
{
    const unsigned nLayerBoundaries (nLayers-1);
    const std::size_t nNodes (layer.getNumberOfNodes());
    const std::vector<MeshLib::Element*> &layer_elements (layer.getElements());
    for (MeshLib::Element* elem : layer_elements)
    {
        const std::size_t nElemNodes (elem->getNumberOfBaseNodes());
        for (unsigned i = 0; i < nElemNodes; ++i)
        {
            if (elem->getNeighbor(i) == nullptr)
            {
                for (unsigned j=0; j<nLayerBoundaries; ++j)
                {
                    const std::size_t offset (j*nNodes);
                    MeshLib::Node* n0 = nodes_[offset + elem->getNodeIndex(i)];
                    MeshLib::Node* n1 = nodes_[offset + elem->getNodeIndex((i+1)%nElemNodes)];
                    MeshLib::Node* n2 = nodes_[offset + nNodes + elem->getNodeIndex((i+1)%nElemNodes)];
                    MeshLib::Node* n3 = nodes_[offset + nNodes + elem->getNodeIndex(i)];

                    if (MathLib::Vector3(*n1, *n2).getLength() > std::numeric_limits<double>::epsilon())
                    {
                        const std::array<MeshLib::Node*,3> tri_nodes = {{ n0, n2, n1 }};
                        elements_.push_back(new MeshLib::Tri(tri_nodes));
                        materials_.push_back(nLayers+j);
                    }
                    if (MathLib::Vector3(*n0, *n3).getLength() > std::numeric_limits<double>::epsilon())
                    {
                        const std::array<MeshLib::Node*,3> tri_nodes = {{ n0, n3, n2 }};
                        elements_.push_back(new MeshLib::Tri(tri_nodes));
                        materials_.push_back(nLayers+j);
                    }
                }
            }
        }
    }
}

void LayeredVolume::removeCongruentElements(std::size_t nLayers, std::size_t nElementsPerLayer)
{
    for (std::size_t i=nLayers-1; i>0; --i)
    {
        const std::size_t lower_offset ((i-1) * nElementsPerLayer);
        const std::size_t upper_offset ( i    * nElementsPerLayer);
        for (std::size_t j=0; j<nElementsPerLayer; ++j)
        {
            MeshLib::Element const*const high (elements_[upper_offset+j]);
            MeshLib::Element *const low  (elements_[lower_offset+j]);

            unsigned count(0);
            const std::size_t nElemNodes (low->getNumberOfBaseNodes());
            for (std::size_t k = 0; k < nElemNodes; ++k)
            {
                if (high->getNodeIndex(k) == low->getNodeIndex(k))
                {
                    low->setNode(k, nodes_[high->getNodeIndex(k)]);
                    count++;
                }
            }

            if (count == nElemNodes)
            {
                delete elements_[upper_offset+j];
                // mark element and material entries for deletion
                elements_[upper_offset+j] = nullptr;
                materials_[upper_offset+j] = -1;
            }
            else
            {
                MeshLib::Node attr = high->getCenterOfGravity();
                attribute_points_.emplace_back(
                    attr[0],
                    attr[1],
                    (attr[2] + low->getCenterOfGravity()[2]) / 2.0,
                    materials_[lower_offset + j]);
            }
        }
    }
    // delete marked entries
    auto elem_vec_end = std::remove(elements_.begin(), elements_.end(), nullptr);
    elements_.erase(elem_vec_end, elements_.end());
    auto mat_vec_end = std::remove(materials_.begin(), materials_.end(), -1);
    materials_.erase(mat_vec_end, materials_.end());
}
