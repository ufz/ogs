/**
 * \file   LayeredVolume.cpp
 * \author Karsten Rink
 * \date   2014-04-11
 * \brief  Implementation of the LayeredVolume class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
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
        return false;

    _elevation_epsilon = calcEpsilon(*rasters[0], *rasters.back());
    if (_elevation_epsilon <= 0)
        return false;

    // remove line elements, only tri + quad remain
    MeshLib::ElementSearch ex(mesh);
    ex.searchByElementType(MeshLib::MeshElemType::LINE);
    std::unique_ptr<MeshLib::Mesh> top(
        removeElements(mesh, ex.getSearchedElementIDs(), "MeshLayer"));
    if (top==nullptr)
        top = std::make_unique<MeshLib::Mesh>(mesh);

    if (!MeshLib::MeshLayerMapper::layerMapping(
            *top, *rasters.back(), noDataReplacementValue))
        return false;

    std::unique_ptr<MeshLib::Mesh> bottom(new MeshLib::Mesh(*top));
    if (!MeshLib::MeshLayerMapper::layerMapping(*bottom, *rasters[0], 0))
        return false;

    this->_minimum_thickness = minimum_thickness;
    _nodes = MeshLib::copyNodeVector(bottom->getNodes());
    _elements = MeshLib::copyElementVector(bottom->getElements(), _nodes);
    if (!_materials.empty())
    {
        ERR("The materials vector is not empty.");
        return false;
    }
    _materials.resize(_elements.size(), 0);

    // map each layer and attach to subsurface mesh
    const std::size_t nRasters (rasters.size());
    for (std::size_t i=1; i<nRasters; ++i)
        this->addLayerToMesh(*top, i, *rasters[i]);

    // close boundaries between layers
    this->addLayerBoundaries(*top, nRasters);
    this->removeCongruentElements(nRasters, top->getNumberOfElements());

    return true;
}

void LayeredVolume::addLayerToMesh(const MeshLib::Mesh &dem_mesh, unsigned layer_id, GeoLib::Raster const& raster)
{
    const std::size_t nNodes (dem_mesh.getNumberOfNodes());
    const std::vector<MeshLib::Node*> &nodes (dem_mesh.getNodes());
    const std::size_t node_id_offset (_nodes.size());
    const std::size_t last_layer_node_offset (node_id_offset-nNodes);

    for (std::size_t i=0; i<nNodes; ++i)
        _nodes.push_back(getNewLayerNode(*nodes[i], *_nodes[last_layer_node_offset + i], raster, _nodes.size()));

    const std::vector<MeshLib::Element*> &layer_elements (dem_mesh.getElements());
    for (MeshLib::Element* elem : layer_elements)
    {
        if (elem->getGeomType() == MeshLib::MeshElemType::TRIANGLE)
        {
            std::array<MeshLib::Node*,3> tri_nodes = {{ _nodes[node_id_offset+elem->getNodeIndex(0)],
                                                        _nodes[node_id_offset+elem->getNodeIndex(1)],
                                                        _nodes[node_id_offset+elem->getNodeIndex(2)] }};
            _elements.push_back(new MeshLib::Tri(tri_nodes));
            _materials.push_back(layer_id);
        }
        else if (elem->getGeomType() == MeshLib::MeshElemType::QUAD)
        {
            std::array<MeshLib::Node*,4> quad_nodes = {{ _nodes[node_id_offset+elem->getNodeIndex(0)],
                                                         _nodes[node_id_offset+elem->getNodeIndex(1)],
                                                         _nodes[node_id_offset+elem->getNodeIndex(2)],
                                                         _nodes[node_id_offset+elem->getNodeIndex(3)] }};
            _elements.push_back(new MeshLib::Quad(quad_nodes));
            _materials.push_back(layer_id);
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
        for (unsigned i=0; i<nElemNodes; ++i)
            if (elem->getNeighbor(i) == nullptr)
                for (unsigned j=0; j<nLayerBoundaries; ++j)
                {
                    const std::size_t offset (j*nNodes);
                    MeshLib::Node* n0 = _nodes[offset + elem->getNodeIndex(i)];
                    MeshLib::Node* n1 = _nodes[offset + elem->getNodeIndex((i+1)%nElemNodes)];
                    MeshLib::Node* n2 = _nodes[offset + nNodes + elem->getNodeIndex((i+1)%nElemNodes)];
                    MeshLib::Node* n3 = _nodes[offset + nNodes + elem->getNodeIndex(i)];

                    if (MathLib::Vector3(*n1, *n2).getLength() > std::numeric_limits<double>::epsilon())
                    {
                        const std::array<MeshLib::Node*,3> tri_nodes = {{ n0, n2, n1 }};
                        _elements.push_back(new MeshLib::Tri(tri_nodes));
                        _materials.push_back(nLayers+j);
                    }
                    if (MathLib::Vector3(*n0, *n3).getLength() > std::numeric_limits<double>::epsilon())
                    {
                        const std::array<MeshLib::Node*,3> tri_nodes = {{ n0, n3, n2 }};
                        _elements.push_back(new MeshLib::Tri(tri_nodes));
                        _materials.push_back(nLayers+j);
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
            MeshLib::Element const*const high (_elements[upper_offset+j]);
            MeshLib::Element *const low  (_elements[lower_offset+j]);

            unsigned count(0);
            const std::size_t nElemNodes (low->getNumberOfBaseNodes());
            for (std::size_t k=0; k<nElemNodes; ++k)
                if (high->getNodeIndex(k) == low->getNodeIndex(k))
                {
                    low->setNode(k, _nodes[high->getNodeIndex(k)]);
                    count++;
                }

            if (count == nElemNodes)
            {
                delete _elements[upper_offset+j];
                // mark element and material entries for deletion
                _elements[upper_offset+j] = nullptr;
                _materials[upper_offset+j] = -1;
            }
            else
            {
                MeshLib::Node attr = high->getCenterOfGravity();
                _attribute_points.emplace_back(
                    attr[0],
                    attr[1],
                    (attr[2] + low->getCenterOfGravity()[2]) / 2.0,
                    _materials[lower_offset + j]);
            }
        }
    }
    // delete marked entries
    auto elem_vec_end = std::remove(_elements.begin(), _elements.end(), nullptr);
    _elements.erase(elem_vec_end, _elements.end());
    auto mat_vec_end = std::remove(_materials.begin(), _materials.end(), -1);
    _materials.erase(mat_vec_end, _materials.end());
}
