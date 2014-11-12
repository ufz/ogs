/**
 * \file   LayeredVolume.cpp
 * \author Karsten Rink
 * \date   2014-04-11
 * \brief  Implementation of the LayeredVolume class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LayeredVolume.h"

#include "Vector3.h"

#include "Raster.h"

#include "MeshLib/Elements/Tri.h"
#include "MeshLib/Elements/Quad.h"
#include "MeshLib/MeshEditing/ElementExtraction.h"
#include "MeshLib/MeshEditing/DuplicateMeshComponents.h"
#include "MeshLib/MeshGenerators/MeshLayerMapper.h"


bool LayeredVolume::createRasterLayers(const MeshLib::Mesh &mesh, const std::vector<GeoLib::Raster const*> &rasters, double noDataReplacementValue)
{
	if (mesh.getDimension() != 2)
		return false;

	_elevation_epsilon = calcEpsilon(*rasters[0], *rasters.back());
	if (_elevation_epsilon <= 0)
		return false;
	
	// remove line elements, only tri + quad remain
	MeshLib::ElementExtraction ex(mesh);
	ex.searchByElementType(MeshElemType::LINE);
	MeshLib::Mesh* top (ex.removeMeshElements("MeshLayer"));
	if (top==nullptr)
		top = new MeshLib::Mesh(mesh);

	if (!MeshLayerMapper::layerMapping(*top, *rasters.back(), noDataReplacementValue))
		return false;

	MeshLib::Mesh* bottom (new MeshLib::Mesh(*top));
	if (!MeshLayerMapper::layerMapping(*bottom, *rasters[0], 0))
	{
		delete top;
		return false;
	}

	_nodes = MeshLib::copyNodeVector(bottom->getNodes());
	_elements = MeshLib::copyElementVector(bottom->getElements(), _nodes);
	delete bottom;

	// map each layer and attach to subsurface mesh
	const std::size_t nRasters (rasters.size());
	for (std::size_t i=1; i<nRasters; ++i)
		this->addLayerToMesh(*top, i, *rasters[i]);

	// close boundaries between layers
	this->addLayerBoundaries(*top, nRasters);
	this->removeCongruentElements(nRasters, top->getNElements());
	delete top;
	return true;
}

void LayeredVolume::addLayerToMesh(const MeshLib::Mesh &dem_mesh, unsigned layer_id, GeoLib::Raster const& raster)
{
	const std::size_t nNodes (dem_mesh.getNNodes());
	const std::vector<MeshLib::Node*> &nodes (dem_mesh.getNodes());
	const std::size_t node_id_offset (_nodes.size());
	const std::size_t last_layer_node_offset (node_id_offset-nNodes);

	for (std::size_t i=0; i<nNodes; ++i)
		_nodes.push_back(getNewLayerNode(*nodes[i], *_nodes[last_layer_node_offset + i], raster, _nodes.size(), _elevation_epsilon));

	const std::vector<MeshLib::Element*> &layer_elements (dem_mesh.getElements());
	for (MeshLib::Element* elem : layer_elements)
	{
		if (elem->getGeomType() == MeshElemType::TRIANGLE)
		{
			std::array<MeshLib::Node*,3> tri_nodes = {{ _nodes[node_id_offset+elem->getNodeIndex(0)],
			                                            _nodes[node_id_offset+elem->getNodeIndex(1)],
			                                            _nodes[node_id_offset+elem->getNodeIndex(2)] }};
			_elements.push_back(new MeshLib::Tri(tri_nodes, layer_id));
		}
		else if (elem->getGeomType() == MeshElemType::QUAD)
		{
			std::array<MeshLib::Node*,4> quad_nodes = {{ _nodes[node_id_offset+elem->getNodeIndex(0)],
			                                             _nodes[node_id_offset+elem->getNodeIndex(1)],
			                                             _nodes[node_id_offset+elem->getNodeIndex(2)],
			                                             _nodes[node_id_offset+elem->getNodeIndex(3)] }};
			_elements.push_back(new MeshLib::Quad(quad_nodes, layer_id));
		}
	}
}

void LayeredVolume::addLayerBoundaries(const MeshLib::Mesh &layer, std::size_t nLayers)
{
	const unsigned nLayerBoundaries (nLayers-1);
	const std::size_t nNodes (layer.getNNodes());
	const std::vector<MeshLib::Element*> &layer_elements (layer.getElements());
	for (MeshLib::Element* elem : layer_elements)
	{
		const std::size_t nElemNodes (elem->getNNodes());
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
						_elements.push_back(new MeshLib::Tri(tri_nodes, nLayers+1+j));
					}
					if (MathLib::Vector3(*n0, *n3).getLength() > std::numeric_limits<double>::epsilon())
					{
						const std::array<MeshLib::Node*,3> tri_nodes = {{ n0, n3, n2 }};
						_elements.push_back(new MeshLib::Tri(tri_nodes, nLayers+1+j));
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
			const std::size_t nElemNodes (low->getNNodes());
			for (std::size_t k=0; k<nElemNodes; ++k)
				if (high->getNodeIndex(k) == low->getNodeIndex(k))
				{
					low->setNode(k, _nodes[high->getNodeIndex(k)]);
					count++;
				}

			if (count == nElemNodes)
			{
				delete _elements[upper_offset+j];
				_elements[upper_offset+j] = nullptr;
			}
			else
			{
				MeshLib::Node attr = high->getCenterOfGravity();
				_attribute_points.push_back(MeshLib::Node(attr[0], attr[1], (attr[2] + low->getCenterOfGravity()[2])/2.0, low->getValue()));
			}
		}
	}
	auto elem_vec_end = std::remove(_elements.begin(), _elements.end(), nullptr);
	_elements.erase(elem_vec_end, _elements.end());
}
