/**
 * \file   LayerVolumes.cpp
 * \author Karsten Rink
 * \date   2014-04-11
 * \brief  Implementation of the LayerVolumes class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LayerVolumes.h"

#include <fstream>
#include <numeric>

#include "GEOObjects.h"
#include "PointVec.h"
#include "Mesh.h"
#include "Node.h"
#include "convertMeshToGeo.h"
#include "Elements/Element.h"
#include "Elements/Tri.h"
#include "Elements/Quad.h"
#include "MeshGenerators/MeshLayerMapper.h"

LayerVolumes::LayerVolumes()
: _invalid_value(-9999), _mesh(nullptr)
{
}

bool LayerVolumes::createGeoVolumes(const MeshLib::Mesh &mesh, const std::vector<std::string> &raster_paths, double noDataReplacementValue)
{
	if (!allRastersExist(raster_paths))
		return false;

	// DEM
	MeshLib::Mesh mesh_layer (mesh);
	if (!MeshLayerMapper::LayerMapping(mesh_layer, raster_paths[0], 0, 0, noDataReplacementValue))
		return false;
	this->addLayerToMesh(mesh_layer, 0);

	// subsurface layers
	const std::size_t nRasters (raster_paths.size());
	for (size_t i=1; i<nRasters; ++i)
	{
		if (!MeshLayerMapper::LayerMapping(mesh_layer, raster_paths[0], 0, 0, _invalid_value))
		{
			this->cleanUpOnError();
			return false;
		}
		this->addLayerToMesh(mesh_layer, i);
	}
	this->addLayerBoundaries(mesh_layer, nRasters);
	_mesh = new MeshLib::Mesh("BoundaryMesh", _nodes, _elements);

	return true;
}

void LayerVolumes::addLayerToMesh(const MeshLib::Mesh &mesh_layer, unsigned layer_id)
{
	const std::vector<MeshLib::Node*> &layer_nodes (mesh_layer.getNodes());
	const std::size_t nNodes (layer_nodes.size());
	const std::size_t node_id_offset (_nodes.size());
	const std::size_t last_layer_offset (node_id_offset-nNodes);
	
	for (std::size_t i=0; i<nNodes; ++i)
	{
		if (layer_id > 0 &&
		   ((*layer_nodes[i])[2] == _invalid_value || 
		   (*layer_nodes[i])[2] < (*_nodes[(layer_id-1)*nNodes+i])[2]))
			_nodes.push_back(new MeshLib::Node(*_nodes[last_layer_offset+i]));
		_nodes.push_back(new MeshLib::Node(*layer_nodes[i]));
	}

	const std::vector<MeshLib::Element*> &layer_elements (mesh_layer.getElements());
	for (MeshLib::Element* elem : layer_elements)
	{
		if (elem->getGeomType() == MeshElemType::TRIANGLE)
		{
			
			std::array<MeshLib::Node*,3> tri_nodes = { _nodes[node_id_offset+elem->getNode(0)->getID()],
			                                           _nodes[node_id_offset+elem->getNode(1)->getID()],
			                                           _nodes[node_id_offset+elem->getNode(2)->getID()] };
			_elements.push_back(new MeshLib::Tri(tri_nodes, layer_id));
		}
		else if (elem->getGeomType() == MeshElemType::QUAD)
		{			
			std::array<MeshLib::Node*,4> quad_nodes = { _nodes[node_id_offset+elem->getNode(0)->getID()],
			                                            _nodes[node_id_offset+elem->getNode(1)->getID()],
			                                            _nodes[node_id_offset+elem->getNode(2)->getID()],
			                                            _nodes[node_id_offset+elem->getNode(3)->getID()] };
			_elements.push_back(new MeshLib::Quad(quad_nodes, layer_id));
		}
	}
}

void LayerVolumes::addLayerBoundaries(const MeshLib::Mesh &layer, std::size_t nLayers)
{
	const unsigned nLayerBoundaries (nLayers-1);
	const std::size_t nNodes (layer.getNNodes());
	const std::vector<MeshLib::Element*> &layer_elements (layer.getElements());
	for (MeshLib::Element* elem : layer_elements)
	{
		for (unsigned i=0; i<3; ++i)
			if (elem->getNeighbor(i) == nullptr)
				for (unsigned j=0; j<nLayerBoundaries; ++j)
				{
					const std::size_t offset (j*nNodes);
					std::array<MeshLib::Node*,4> quad_nodes = 
						{ 
						  _nodes[offset + elem->getNode(i)->getID()],
						  _nodes[offset + elem->getNode((i+1)%3)->getID()],
						  _nodes[offset + nNodes + elem->getNode((i+1)%3)->getID()],
						  _nodes[offset + nNodes + elem->getNode(i)->getID()] 
						};
					_elements.push_back(new MeshLib::Quad(quad_nodes, nLayers+j));
				}
	}
}

bool LayerVolumes::addGeometry(GeoLib::GEOObjects &geo_objects) const
{
	if (_mesh == nullptr)
		return false;
	MeshLib::convertMeshToGeo(*_mesh, geo_objects);
	return true;
}

bool LayerVolumes::allRastersExist(const std::vector<std::string> &raster_paths) const
{
	for (auto raster = raster_paths.begin(); raster != raster_paths.end(); ++raster)
	{
		std::ifstream file_stream (*raster, std::ifstream::in);
		if (!file_stream.good())
			return false;
		file_stream.close();
	}
	return true;
}

void LayerVolumes::cleanUpOnError()
{
	std::for_each(_nodes.begin(), _nodes.end(), [](MeshLib::Node *node) { delete node; });
	std::for_each(_elements.begin(), _elements.end(), [](MeshLib::Element *elem) { delete elem; });
}
