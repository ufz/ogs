/**
 * \file   MeshLayerMapper.cpp
 * \author Karsten Rink
 * \date   2010-11-01
 * \brief  Implementation of the MeshLayerMapper class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshLayerMapper.h"

#include <algorithm>

#include "logog/include/logog.hpp"

#include "FileIO/AsciiRasterInterface.h"

#include "Raster.h"

#include "Elements/Tet.h"
#include "Elements/Hex.h"
#include "Elements/Pyramid.h"
#include "Elements/Prism.h"
#include "MeshSurfaceExtraction.h"

#include "MathTools.h"

const unsigned MeshLayerMapper::_pyramid_base[3][4] =
{
	{1, 3, 4, 2}, // Point 4 missing
	{2, 4, 3, 0}, // Point 5 missing
	{0, 3, 4, 1}, // Point 6 missing
};

MeshLib::Mesh* MeshLayerMapper::createStaticLayers(MeshLib::Mesh const& mesh, std::vector<float> const& layer_thickness_vector, std::string const& mesh_name) const
{
	std::vector<float> thickness;
	for (std::size_t i=0; i<layer_thickness_vector.size(); ++i)
		if (layer_thickness_vector[i] > std::numeric_limits<float>::epsilon())
			thickness.push_back(layer_thickness_vector[i]);
		else
			WARN ("Ignoring layer %d with thickness %f.", i, layer_thickness_vector[i]);

	const std::size_t nLayers(thickness.size());
	if (nLayers < 1 || mesh.getDimension() != 2)
	{
		ERR("MeshLayerMapper::createStaticLayers(): A 2D mesh with nLayers > 0 is required as input.");
		return nullptr;
	}

	const std::size_t nNodes = mesh.getNNodes();
	// count number of 2d elements in the original mesh
	const std::size_t nElems (std::count_if(mesh.getElements().begin(), mesh.getElements().end(),
			[](MeshLib::Element const* elem) { return (elem->getDimension() == 2);}));

	const std::size_t nOrgElems (mesh.getNElements());
	const std::vector<MeshLib::Node*> &nodes = mesh.getNodes();
	const std::vector<MeshLib::Element*> &elems = mesh.getElements();
	std::vector<MeshLib::Node*> new_nodes(nNodes + (nLayers * nNodes));
	std::vector<MeshLib::Element*> new_elems;
	new_elems.reserve(nElems * nLayers);
	double z_offset (0.0);

	for (unsigned layer_id = 0; layer_id <= nLayers; ++layer_id)
	{
		// add nodes for new layer
		unsigned node_offset (nNodes * layer_id);
		if (layer_id > 0) z_offset += thickness[layer_id-1];

		std::transform(nodes.cbegin(), nodes.cend(), new_nodes.begin() + node_offset,
			[&z_offset](MeshLib::Node* node){ return new MeshLib::Node((*node)[0], (*node)[1], (*node)[2]-z_offset); });

		// starting with 2nd layer create prism or hex elements connecting the last layer with the current one
		if (layer_id == 0)
			continue;

		node_offset -= nNodes;
		const unsigned mat_id (nLayers - layer_id);

		for (unsigned i = 0; i < nOrgElems; ++i)
		{
			const MeshLib::Element* sfc_elem( elems[i] );
			if (sfc_elem->getDimension() < 2) // ignore line-elements
				continue;
				
			const unsigned nElemNodes(sfc_elem->getNNodes());
			MeshLib::Node** e_nodes = new MeshLib::Node*[2*nElemNodes];

			for (unsigned j=0; j<nElemNodes; ++j)
			{
				const unsigned node_id = sfc_elem->getNode(j)->getID() + node_offset;
				e_nodes[j] = new_nodes[node_id+nNodes];
				e_nodes[j+nElemNodes] = new_nodes[node_id];
			}
			if (sfc_elem->getGeomType() == MeshElemType::TRIANGLE)	// extrude triangles to prism
				new_elems.push_back (new MeshLib::Prism(e_nodes, mat_id));
			else if (sfc_elem->getGeomType() == MeshElemType::QUAD)	// extrude quads to hexes
				new_elems.push_back (new MeshLib::Hex(e_nodes, mat_id));
		}
	}
	return new MeshLib::Mesh(mesh_name, new_nodes, new_elems);
}

bool MeshLayerMapper::createRasterLayers(MeshLib::Mesh const& mesh, std::vector<GeoLib::Raster const*> const& rasters, double noDataReplacementValue)
{
	const std::size_t nLayers(rasters.size());
	if (nLayers < 2 || mesh.getDimension() != 2)
	{
		ERR("MeshLayerMapper::createRasterLayers(): A 2D mesh and at least two rasters required as input.");
		return false;
	}

	MeshLib::Mesh* top (new MeshLib::Mesh(mesh));
	if (!layerMapping(*top, *rasters.back(), 0))
		return false;

	MeshLib::Mesh* bottom (new MeshLib::Mesh(mesh));
	if (!layerMapping(*bottom, *rasters[0], 0))
	{
		delete top;
		return false;
	}

	std::size_t const nNodes = mesh.getNNodes();
	_nodes.reserve(nLayers * nNodes);

	// number of triangles in the original mesh
	std::size_t const nElems (std::count_if(mesh.getElements().begin(), mesh.getElements().end(),
		[](MeshLib::Element const* elem) { return (elem->getGeomType() == MeshElemType::TRIANGLE);}));
	_elements.reserve(nElems * (nLayers-1));

	// add bottom layer
	std::vector<MeshLib::Node*> const& nodes = bottom->getNodes();
	for (MeshLib::Node* node : nodes)
		_nodes.push_back(new MeshLib::Node(*node));
	delete bottom;
	
	// add the other layers
	for (std::size_t i=1; i<nLayers; ++i)
		addLayerToMesh(*top, i, *rasters[i]);

	delete top;
	return true;
}

void MeshLayerMapper::addLayerToMesh(const MeshLib::Mesh &dem_mesh, unsigned layer_id, GeoLib::Raster const& raster)
{
    std::size_t const nNodes = dem_mesh.getNNodes();
    std::vector<MeshLib::Node*> const& nodes = dem_mesh.getNodes();
    int const last_layer_node_offset = (layer_id-1) * nNodes;

    // add nodes for new layer
    for (std::size_t i=0; i<nNodes; ++i)
        _nodes.push_back(getNewLayerNode(*nodes[i], *_nodes[last_layer_node_offset + i], raster, _nodes.size(), _elevation_epsilon));

    std::vector<MeshLib::Element*> const& elems = dem_mesh.getElements();
    std::size_t const nElems (dem_mesh.getNElements());

    for (std::size_t i=0; i<nElems; ++i)
    {
        MeshLib::Element* elem (elems[i]);
        if (elem->getGeomType() != MeshElemType::TRIANGLE)
            continue;
        unsigned node_counter(3), missing_idx(0);
        std::array<MeshLib::Node*, 6> new_elem_nodes;
        for (unsigned j=0; j<3; ++j)
        {
            new_elem_nodes[j] = _nodes[_nodes[last_layer_node_offset + elem->getNodeIndex(j)]->getID()];
            new_elem_nodes[node_counter] = (_nodes[last_layer_node_offset + elem->getNodeIndex(j) + nNodes]);
            if (new_elem_nodes[j]->getID() != new_elem_nodes[node_counter]->getID())
                node_counter++;
            else
                missing_idx = j;
        }

        switch (node_counter)
        {
        case 6:
            _elements.push_back(new MeshLib::Prism(new_elem_nodes, layer_id-1));
            break;
        case 5:
            std::array<MeshLib::Node*, 5> pyramid_nodes;
            pyramid_nodes[0] = new_elem_nodes[_pyramid_base[missing_idx][0]];
            pyramid_nodes[1] = new_elem_nodes[_pyramid_base[missing_idx][1]];
            pyramid_nodes[2] = new_elem_nodes[_pyramid_base[missing_idx][2]];
            pyramid_nodes[3] = new_elem_nodes[_pyramid_base[missing_idx][3]];
            pyramid_nodes[4] = new_elem_nodes[missing_idx];
            _elements.push_back(new MeshLib::Pyramid(pyramid_nodes, layer_id-1));
            break;
        case 4:
            std::array<MeshLib::Node*, 4> tet_nodes;
            std::copy(new_elem_nodes.begin(), new_elem_nodes.begin() + node_counter, tet_nodes.begin());
            _elements.push_back(new MeshLib::Tet(tet_nodes, layer_id-1));
            break;
        default:
            continue;
        }
    }
}

bool MeshLayerMapper::layerMapping(MeshLib::Mesh &new_mesh, std::string const& rasterfile, double noDataReplacementValue = 0.0)
{
	const GeoLib::Raster *raster(FileIO::AsciiRasterInterface::getRasterFromASCFile(rasterfile));
	if (! raster) {
		ERR("MshLayerMapper::layerMapping - could not read raster file %s", rasterfile.c_str());
		return false;
	}
	const bool result = layerMapping(new_mesh, *raster, noDataReplacementValue);
	delete raster;
	return result;
}

bool MeshLayerMapper::layerMapping(MeshLib::Mesh &new_mesh, GeoLib::Raster const& raster, double noDataReplacementValue = 0.0)
{
	if (new_mesh.getDimension() != 2)
	{
		ERR("MshLayerMapper::layerMapping - requires 2D mesh");
		return false;
	}

	const double x0(raster.getOrigin()[0]);
	const double y0(raster.getOrigin()[1]);
	const double delta(raster.getRasterPixelSize());
	const double no_data(raster.getNoDataValue());

	const std::pair<double, double> xDim(x0, x0 + raster.getNCols() * delta); // extension in x-dimension
	const std::pair<double, double> yDim(y0, y0 + raster.getNRows() * delta); // extension in y-dimension

	const size_t nNodes (new_mesh.getNNodes());
	const std::vector<MeshLib::Node*> &nodes = new_mesh.getNodes();
	for (unsigned i = 0; i < nNodes; ++i)
	{
		if (!raster.isPntOnRaster(*nodes[i]))
		{
			// use either default value or elevation from layer above
			nodes[i]->updateCoordinates((*nodes[i])[0], (*nodes[i])[1], noDataReplacementValue);
			continue;
		}

		double elevation (raster.interpolateValueAtPoint(*nodes[i]));
		if (std::abs(elevation - no_data) < std::numeric_limits<double>::epsilon()) 
			elevation = noDataReplacementValue;
		nodes[i]->updateCoordinates((*nodes[i])[0], (*nodes[i])[1], elevation);
	}

	return true;
}
