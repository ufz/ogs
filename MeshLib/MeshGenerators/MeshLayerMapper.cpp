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

// stl
#include <algorithm>

// ThirdParty/logog
#include "logog/include/logog.hpp"

#include "MeshLayerMapper.h"

#include "FileIO/AsciiRasterInterface.h"

// GeoLib
#include "Raster.h"

#include "Mesh.h"
#include "Node.h"
#include "Elements/Element.h"
#include "Elements/Tet.h"
#include "Elements/Hex.h"
#include "Elements/Pyramid.h"
#include "Elements/Prism.h"
#include "MeshSurfaceExtraction.h"
#include "MeshQuality/MeshValidation.h"

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

MeshLib::Mesh* MeshLayerMapper::createRasterLayers(MeshLib::Mesh const& mesh, std::vector<std::string> const& raster_paths, std::string const& mesh_name) const
{
	if (mesh.getDimension() != 2 || !allRastersExist(raster_paths))
		return false;

	std::vector<GeoLib::Raster const*> rasters;
	rasters.reserve(raster_paths.size());
	for (auto path = raster_paths.begin(); path != raster_paths.end(); ++path)
		rasters.push_back(FileIO::AsciiRasterInterface::getRasterFromASCFile(*path));

	MeshLib::Mesh* result = this->createRasterLayers(mesh, rasters, mesh_name);
	std::for_each(rasters.begin(), rasters.end(), [](GeoLib::Raster const*const raster){ delete raster; });
	return result;
}

MeshLib::Mesh* MeshLayerMapper::createRasterLayers(MeshLib::Mesh const& mesh, std::vector<GeoLib::Raster const*> const& rasters, std::string const& mesh_name) const
{
    const std::size_t nLayers(rasters.size());
	if (nLayers < 1 || mesh.getDimension() != 2)
	{
		ERR("MeshLayerMapper::createStaticLayers(): A 2D mesh with nLayers > 0 is required as input.");
		return nullptr;
	}

    MeshLib::Mesh* dem_mesh (new MeshLib::Mesh(mesh));
    if (layerMapping(*dem_mesh, *rasters.back(), 0, 0, 0))
    {
        std::size_t const nNodes = mesh.getNNodes();
	    std::vector<MeshLib::Node*> new_nodes;
        new_nodes.reserve(nLayers * nNodes);

	    // number of triangles in the original mesh
	    std::size_t const nElems (std::count_if(mesh.getElements().begin(), mesh.getElements().end(),
			[](MeshLib::Element const* elem) { return (elem->getGeomType() == MeshElemType::TRIANGLE);}));
	    std::vector<MeshLib::Element*> new_elems;
        new_elems.reserve(nElems * (nLayers-1));

        for (std::size_t i=0; i<nLayers; ++i)
            addLayerToMesh(*dem_mesh, i, *rasters[i], new_nodes, new_elems);

        MeshLib::Mesh* result = new MeshLib::Mesh(mesh_name, new_nodes, new_elems);
        MeshLib::MeshValidation::removeUnusedMeshNodes(*result);
        return result;

    }
}

void MeshLayerMapper::addLayerToMesh(const MeshLib::Mesh &dem_mesh, unsigned layer_id, GeoLib::Raster const& raster, std::vector<MeshLib::Node*> &new_nodes, std::vector<MeshLib::Element*> &new_elems) const
{
    std::size_t const nNodes = dem_mesh.getNNodes();
	std::vector<MeshLib::Node*> const& nodes = dem_mesh.getNodes();
    int const last_layer_node_offset = (layer_id-1) * nNodes;
    unsigned const layer_node_offset = layer_id * nNodes;

	// add nodes for new layer
    for (std::size_t i=0; i<nNodes; ++i)
    {
        // min of dem elevation and layer elevation
        double const elevation = std::min(raster.interpolateValueAtPoint(*nodes[i]), (*nodes[i])[2]);

        if ((layer_id == 0) || 
            (elevation - (*new_nodes[last_layer_node_offset + i])[2] > std::numeric_limits<double>::epsilon()))
            new_nodes.push_back(new MeshLib::Node((*nodes[i])[0], (*nodes[i])[1], elevation, (layer_id * nNodes) + i));
        else
            new_nodes.push_back(new MeshLib::Node((*nodes[i])[0], (*nodes[i])[1], elevation, new_nodes[last_layer_node_offset +i]->getID()));
    }

    if (layer_id == 0)
        return;

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
            new_elem_nodes[j] = new_nodes[last_layer_node_offset + elem->getNodeIndex(j)];
            new_elem_nodes[node_counter] = (new_nodes[last_layer_node_offset + elem->getNodeIndex(j) + nNodes]);
            if (new_elem_nodes[j]->getID() != new_elem_nodes[node_counter]->getID())
                node_counter++;
            else
                missing_idx = j;
        }

        switch (node_counter)
        {
        case 6:
            new_elems.push_back(new MeshLib::Prism(new_elem_nodes, layer_id));
            break;
        case 5:
            std::array<MeshLib::Node*, 5> pyramid_nodes;
            pyramid_nodes[0] = new_elem_nodes[_pyramid_base[missing_idx][0]];
            pyramid_nodes[1] = new_elem_nodes[_pyramid_base[missing_idx][1]];
            pyramid_nodes[2] = new_elem_nodes[_pyramid_base[missing_idx][2]];
            pyramid_nodes[3] = new_elem_nodes[_pyramid_base[missing_idx][3]];
            pyramid_nodes[4] = new_elem_nodes[missing_idx];
            new_elems.push_back(new MeshLib::Pyramid(pyramid_nodes, layer_id));
            break;
        case 4:
            std::array<MeshLib::Node*, 4> tet_nodes;
            std::copy(new_elem_nodes.begin(), new_elem_nodes.begin() + node_counter, tet_nodes.begin());
            new_elems.push_back(new MeshLib::Tet(tet_nodes, layer_id));
            break;
        default:
            continue;
        }
    }
}


bool MeshLayerMapper::layerMapping(MeshLib::Mesh &new_mesh, const std::string &rasterfile,
                                   const unsigned nLayers, const unsigned layer_id, double noDataReplacementValue = 0.0) const
{
	const GeoLib::Raster *raster(FileIO::AsciiRasterInterface::getRasterFromASCFile(rasterfile));
	if (! raster) {
		ERR("MshLayerMapper::LayerMapping - could not read raster file %s", rasterfile.c_str());
		return false;
	}
	const bool result = layerMapping(new_mesh, *raster, nLayers, layer_id, noDataReplacementValue);
	delete raster;
	return result;
}

bool MeshLayerMapper::layerMapping(MeshLib::Mesh &new_mesh, const GeoLib::Raster &raster,
                                   const unsigned nLayers, const unsigned layer_id, double noDataReplacementValue = 0.0) const
{
	if (nLayers < layer_id)
	{
		ERR("MshLayerMapper::LayerMapping() - Mesh has only %d Layers, cannot assign layer %d.", nLayers, layer_id);
		return false;
	}

	const double x0(raster.getOrigin()[0]);
	const double y0(raster.getOrigin()[1]);
	const double delta(raster.getRasterPixelSize());
	const double no_data(raster.getNoDataValue());
	double const*const elevation(raster.begin());

	const std::pair<double, double> xDim(x0, x0 + raster.getNCols() * delta); // extension in x-dimension
	const std::pair<double, double> yDim(y0, y0 + raster.getNRows() * delta); // extension in y-dimension

	const size_t nNodes (new_mesh.getNNodes());
	const size_t nNodesPerLayer (nNodes / (nLayers+1));

	const size_t firstNode (layer_id * nNodesPerLayer);
	const size_t lastNode  (firstNode + nNodesPerLayer);

	const double half_delta = 0.5*delta;
	const std::vector<MeshLib::Node*> &nodes = new_mesh.getNodes();
	for (unsigned i = firstNode; i < lastNode; ++i)
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

MeshLib::Mesh* MeshLayerMapper::blendLayersWithSurface(MeshLib::Mesh &mesh, const unsigned nLayers, const std::string &dem_raster) const
{
	// construct surface mesh from DEM
	const MathLib::Vector3 dir(0,0,1);
	MeshLib::Mesh* dem = MeshLib::MeshSurfaceExtraction::getMeshSurface(mesh, dir);
	this->layerMapping(*dem, dem_raster, 0, 0);
	const std::vector<MeshLib::Node*> &dem_nodes (dem->getNodes());

	const std::vector<MeshLib::Node*> &mdl_nodes (mesh.getNodes());
	const size_t nNodes (mesh.getNNodes());
	const size_t nNodesPerLayer = nNodes / (nLayers+1);
	std::vector<bool> is_surface_node(nNodes, false);
	std::vector<bool> nodes_below_surface(nNodes, false);

	// check if bottom layer nodes are below DEM
	const unsigned bottom_firstNode = nLayers * nNodesPerLayer;
	const unsigned bottom_lastNode  = bottom_firstNode + nNodesPerLayer;
	for (unsigned i = bottom_firstNode; i < bottom_lastNode; ++i)
	{
		nodes_below_surface[i]=true;
		const double* coords = mdl_nodes[i]->getCoords();
		const double* dem_coords = dem_nodes[i-bottom_firstNode]->getCoords();
		if (coords[2] >= dem_coords[2])
		{
			WARN("Node %d (in bottom-layer) is above surface node %d. (%f, %f)", i, (i-bottom_firstNode), coords[2], dem_coords[2]);
			is_surface_node[i] = true;
		}
	}

	// for all other layers:
	// if node < dem-node: do nothing
	// if node > dem-node:
	//		if first node above surface: map to dem and mark as surface node
	//		else remove node (i.e. don't copy it)
	for (int layer_id=nLayers-1; layer_id>=0; --layer_id)
	{
		const size_t firstNode = layer_id * nNodesPerLayer;
		const size_t lastNode  = firstNode + nNodesPerLayer;

		for(unsigned i = firstNode; i < lastNode; ++i)
		{
			if (is_surface_node[i+nNodesPerLayer])
				is_surface_node[i]=true;
			else
			{
				nodes_below_surface[i]=true;
				MeshLib::Node* node (mdl_nodes[i]);
				const double* coords = node->getCoords();
				const double* dem_coords = dem_nodes[i-firstNode]->getCoords();
				if (coords[2] > dem_coords[2])
				{
					node->updateCoordinates(dem_coords[0], dem_coords[1], dem_coords[2]);
					is_surface_node[i] = true;
				}
			}
		}
	}
	delete dem; // no longer needed

	// copy valid nodes to new node vector
	std::vector<MeshLib::Node*> new_nodes;
	std::vector<int> node_index_map(nNodes, -1);
	size_t node_count(0);
	for (unsigned j=0; j<nNodes; ++j)
		if (nodes_below_surface[j])
		{
			new_nodes.push_back(new MeshLib::Node(mdl_nodes[j]->getCoords(), mdl_nodes[j]->getID()));
			node_index_map[j]=node_count++;
		}

	// copy elements (old elements need to have at least 4 nodes remaining and form a 3d element
	const std::vector<MeshLib::Element*> &mdl_elements (mesh.getElements());
	const size_t nElems (mesh.getNElements());
	std::vector<MeshLib::Element*> new_elements;
	for (unsigned j=0; j<nElems; ++j)
	{
		const MeshLib::Element* elem = mdl_elements[j];

		size_t count(0);
		for (unsigned i=0; i<6; ++i) // check top surface of prism
			if (nodes_below_surface[elem->getNode(i)->getID()]) ++count;

		if (count==6) // copy prism elements if all six nodes are valid
		{
			MeshLib::Node** e_nodes = new MeshLib::Node*[count];
			for (unsigned i=0; i<6; ++i)
				e_nodes[i] = new_nodes[node_index_map[elem->getNode(i)->getID()]];

			MeshLib::Element* prism (new MeshLib::Prism(e_nodes, elem->getValue()));
			new_elements.push_back(prism);
		}
		else if (count==5) // change the current element to two tetrahedra if only five nodes are valid
		{
			MeshLib::Node** e_nodes = new MeshLib::Node*[count];
			unsigned top_idx(6);
			for (unsigned i=3; i<6; ++i) // find node that has been cut
				if (!nodes_below_surface[elem->getNode(i)->getID()])
					top_idx = i-3;

			// construct pyramid element based on missing node
			unsigned idx1 ((top_idx+1)%3);
			unsigned idx2 ((top_idx+2)%3);
			e_nodes[0] = new_nodes[node_index_map[elem->getNode(idx1)->getID()]];
			e_nodes[1] = new_nodes[node_index_map[elem->getNode(idx1+3)->getID()]];
			e_nodes[2] = new_nodes[node_index_map[elem->getNode(idx2+3)->getID()]];
			e_nodes[3] = new_nodes[node_index_map[elem->getNode(idx2)->getID()]];
			e_nodes[4] = new_nodes[node_index_map[elem->getNode(top_idx)->getID()]];

			MeshLib::Element* pyr (new MeshLib::Pyramid(e_nodes, elem->getValue()));
			new_elements.push_back(pyr);
		}
		else if (count==4) // change the current element to a tetrahedron if only four nodes are valid
		{
			MeshLib::Node** e_nodes = new MeshLib::Node*[count];
			for (unsigned i=0; i<3; ++i) // first three nodes are the bottom-face
			{
				unsigned idx (elem->getNode(i)->getID());
				if (nodes_below_surface[idx])
					e_nodes[i] = new_nodes[node_index_map[idx]];
				else
					e_nodes[i] = NULL;
			}

			if (e_nodes[0] && e_nodes[1] && e_nodes[2]) //make sure that the 4 remaining nodes don't form a quad
			{
				for (unsigned i=3; i<6; ++i) // last node
				{
					unsigned idx (elem->getNode(i)->getID());
					if (nodes_below_surface[idx])
					{
						e_nodes[3] = new_nodes[node_index_map[idx]];
						break;
					}
				}

				MeshLib::Element* tet (new MeshLib::Tet(e_nodes, elem->getValue()));
				new_elements.push_back(tet);
			}
			else delete e_nodes;
		}
		// else remove element, if less than four nodes are valid
	}
	return new MeshLib::Mesh("SubsurfaceMesh", new_nodes, new_elements);
}

bool MeshLayerMapper::allRastersExist(const std::vector<std::string> &raster_paths) const
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

