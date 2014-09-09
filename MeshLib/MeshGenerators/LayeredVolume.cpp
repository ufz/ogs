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

#include <fstream>
#include <numeric>

#include "Vector3.h"

#include "GEOObjects.h"
#include "PointVec.h"
#include "Mesh.h"
#include "convertMeshToGeo.h"
#include "Elements/Element.h"
#include "Elements/Tri.h"
#include "Elements/Quad.h"
#include "MeshGenerators/MeshLayerMapper.h"
#include "MeshQuality/MeshValidation.h"
#include "MeshEditing/ElementExtraction.h"


const double LayeredVolume::_invalid_value = -9999;

LayeredVolume::LayeredVolume()
: _elevation_epsilon(0.0001), _mesh(nullptr)
{
}

bool LayeredVolume::createGeoVolumes(const MeshLib::Mesh &mesh, const std::vector<std::string> &raster_paths, double noDataReplacementValue)
{
	if (mesh.getDimension() != 2 || !allRastersExist(raster_paths))
		return false;

	std::vector<GeoLib::Raster const*> rasters;
	rasters.reserve(raster_paths.size());
	for (auto path = raster_paths.begin(); path != raster_paths.end(); ++path)
		rasters.push_back(GeoLib::Raster::getRasterFromASCFile(*path));

	const bool result = this->createGeoVolumes(mesh, rasters, noDataReplacementValue);
	std::for_each(rasters.begin(), rasters.end(), [](GeoLib::Raster const*const raster){ delete raster; });
	return result;
}


bool LayeredVolume::createGeoVolumes(const MeshLib::Mesh &mesh, const std::vector<GeoLib::Raster const*> &rasters, double noDataReplacementValue)
{
	if (mesh.getDimension() != 2)
		return false;

	_elevation_epsilon = calcEpsilon(*rasters[0], *rasters.back());
	if (_elevation_epsilon <= 0)
		return false;

	// remove line elements, only tri + quad remain
	MeshLib::ElementExtraction ex(mesh);
	ex.searchByElementType(MeshElemType::LINE);
	MeshLib::Mesh* mesh_layer (ex.removeMeshElements("MeshLayer"));
	if (mesh_layer==nullptr)
		mesh_layer = new MeshLib::Mesh(mesh);

	// map each layer and attach to subsurface mesh
    MeshLayerMapper const mapper;
	const std::size_t nRasters (rasters.size());
	for (size_t i=0; i<nRasters; ++i)
	{
		const double replacement_value = (i==(nRasters-1)) ? noDataReplacementValue : _invalid_value;
		if (!mapper.layerMapping(*mesh_layer, *rasters[i], 0, 0, replacement_value))
		{
			this->cleanUpOnError();
			return false;
		}
		this->addLayerToMesh(*mesh_layer, i);
	}
	// close boundaries between layers
	this->addLayerBoundaries(*mesh_layer, nRasters);
	this->removeCongruentElements(nRasters, mesh_layer->getNElements());
	delete mesh_layer;
	_mesh = new MeshLib::Mesh("BoundaryMesh", _nodes, _elements);
	MeshLib::MeshValidation::removeUnusedMeshNodes(*_mesh);
	return true;
}

void LayeredVolume::addLayerToMesh(const MeshLib::Mesh &mesh_layer, unsigned layer_id)
{
	const std::vector<MeshLib::Node*> &layer_nodes (mesh_layer.getNodes());
	const std::size_t nNodes (layer_nodes.size());
	const std::size_t node_id_offset (_nodes.size());
	const std::size_t last_layer_offset (node_id_offset-nNodes);

	for (std::size_t i=0; i<nNodes; ++i)
	{
		if (layer_id > 0 &&
		   ((*layer_nodes[i])[2] == _invalid_value ||
		    (*_nodes[last_layer_offset+i])[2]-(*layer_nodes[i])[2] < _elevation_epsilon))
			_nodes.push_back(new MeshLib::Node(*_nodes[last_layer_offset+i]));
		else
			_nodes.push_back(new MeshLib::Node(layer_nodes[i]->getCoords(), _nodes.size()));
	}

	const std::vector<MeshLib::Element*> &layer_elements (mesh_layer.getElements());
	for (MeshLib::Element* elem : layer_elements)
	{
		if (elem->getGeomType() == MeshElemType::TRIANGLE)
		{
			std::array<MeshLib::Node*,3> tri_nodes = {{ _nodes[node_id_offset+elem->getNodeIndex(0)],
			                                           _nodes[node_id_offset+elem->getNodeIndex(1)],
			                                           _nodes[node_id_offset+elem->getNodeIndex(2)] }};
			_elements.push_back(new MeshLib::Tri(tri_nodes, layer_id+1));
		}
		else if (elem->getGeomType() == MeshElemType::QUAD)
		{
			std::array<MeshLib::Node*,4> quad_nodes = {{ _nodes[node_id_offset+elem->getNodeIndex(0)],
			                                            _nodes[node_id_offset+elem->getNodeIndex(1)],
			                                            _nodes[node_id_offset+elem->getNodeIndex(2)],
			                                            _nodes[node_id_offset+elem->getNodeIndex(3)] }};
			_elements.push_back(new MeshLib::Quad(quad_nodes, layer_id+1));
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
	for (std::size_t i=1; i<nLayers; ++i)
	{
		const std::size_t upper_offset ((i-1) * nElementsPerLayer);
		const std::size_t lower_offset ( i    * nElementsPerLayer);
		for (std::size_t j=0; j<nElementsPerLayer; ++j)
		{
			MeshLib::Element const*const high (_elements[upper_offset+j]);
			MeshLib::Element *low  (_elements[lower_offset+j]);

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

bool LayeredVolume::exportToGeometry(GeoLib::GEOObjects &geo_objects) const
{
	if (_mesh == nullptr)
		return false;
	MeshLib::convertMeshToGeo(*_mesh, geo_objects, std::numeric_limits<double>::min());
	return true;
}

double LayeredVolume::calcEpsilon(const GeoLib::Raster &high, const GeoLib::Raster &low)
{
	const double max (*std::max_element(high.begin(), high.end()));
	const double min (*std::min_element( low.begin(),  low.end()));
	return ((max-min)*1e-07);
}

bool LayeredVolume::allRastersExist(const std::vector<std::string> &raster_paths) const
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

void LayeredVolume::cleanUpOnError()
{
	std::for_each(_nodes.begin(), _nodes.end(), [](MeshLib::Node *node) { delete node; });
	std::for_each(_elements.begin(), _elements.end(), [](MeshLib::Element *elem) { delete elem; });
}
