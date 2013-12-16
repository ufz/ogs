/**
 * \file
 * \author Karsten Rink
 * \date   2010-11-01
 * \brief  Implementation of the MshLayerMapper class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// stl
#include <algorithm>

// ThirdParty/logog
#include "logog/include/logog.hpp"

#include "MshLayerMapper.h"
// GeoLib
#include "Raster.h"

#include "Mesh.h"
#include "Node.h"
#include "Elements/Element.h"
#include "Elements/Tet.h"
#include "Elements/Hex.h"
#include "Elements/Pyramid.h"
#include "Elements/Prism.h"
#include "MeshEditing/removeMeshNodes.h"
#include "MeshSurfaceExtraction.h"
#include "MathTools.h"


MeshLib::Mesh* MshLayerMapper::CreateLayers(const MeshLib::Mesh* mesh, const std::vector<float> &thickness)
{
	std::size_t nLayers(thickness.size());
	bool throw_error(false);
	for (unsigned i=0; i<nLayers; ++i)
		if (thickness[i]<=0)
			throw_error = true;
	if (nLayers < 1 || mesh->getDimension() != 2)
		throw_error = true;

	if (throw_error)
	{
		ERR("MshLayerMapper::CreateLayers(): A 2D mesh with nLayers > 0 is required as input.");
		return nullptr;
	}

	const size_t nNodes = mesh->getNNodes();
	// count number of 2d elements in the original mesh
	const std::size_t nElems (std::count_if(mesh->getElements().begin(), mesh->getElements().end(),
			[](MeshLib::Element const* elem) { return (elem->getDimension() == 2);}));

	const std::vector<MeshLib::Node*> nodes = mesh->getNodes();
	const std::vector<MeshLib::Element*> elems = mesh->getElements();
	std::vector<MeshLib::Node*> new_nodes(nNodes + (nLayers * nNodes));
	std::vector<MeshLib::Element*> new_elems(nElems * nLayers);
	double z_offset(0.0);

	for (unsigned layer_id = 0; layer_id <= nLayers; ++layer_id)
	{
		// add nodes for new layer
		unsigned node_offset (nNodes * layer_id);
		unsigned elem_offset (nElems * (layer_id-1));
		if (layer_id>0) z_offset += thickness[layer_id-1];
		for (unsigned i = 0; i < nNodes; ++i)
		{
			const double* coords = nodes[i]->getCoords();
			new_nodes[node_offset+i] = new MeshLib::Node(coords[0], coords[1], coords[2]-z_offset, node_offset+i);
		}

		// starting with 2nd layer create prism or hex elements connecting the last layer with the current one
		if (layer_id > 0)
		{
			if (thickness[layer_id-1] > 0)
			{
				node_offset -= nNodes;
				const unsigned mat_id (nLayers - layer_id);

				// counts the 2d elements (within the layer),
				// used as a part of index computation in new_elems
				std::size_t cnt(0);
				for (unsigned i = 0; i < mesh->getNElements(); ++i)
				{
					const MeshLib::Element* sfc_elem( elems[i] );
					if (sfc_elem->getDimension() == 2)
					{
						const unsigned nElemNodes(sfc_elem->getNNodes());
						MeshLib::Node** e_nodes = new MeshLib::Node*[2*nElemNodes];

						for (unsigned j=0; j<nElemNodes; ++j)
						{
							const unsigned node_id = sfc_elem->getNode(j)->getID() + node_offset;
							e_nodes[j] = new_nodes[node_id+nNodes];
							e_nodes[j+nElemNodes] = new_nodes[node_id];
						}
						if (sfc_elem->getGeomType() == MeshElemType::TRIANGLE)	// extrude triangles to prism
							new_elems[elem_offset+cnt] = new MeshLib::Prism(e_nodes, mat_id);
						else if (sfc_elem->getGeomType() == MeshElemType::QUAD)	// extrude quads to hexes
							new_elems[elem_offset+cnt] = new MeshLib::Hex(e_nodes, mat_id);
						cnt++;
					}
					else
					{
						WARN("MshLayerMapper::CreateLayers() - Method can only handle 2D mesh elements.");
						WARN("Skipping Element %d of type \"%s\".", i, MeshElemType2String(sfc_elem->getGeomType()).c_str());
					}
				}
			}
			else
			{
				ERR("Error in MshLayerMapper::CreateLayers() - Layer thickness for layer %d is %f (needs to be >0).", (layer_id-1), thickness[layer_id-1]);
				return nullptr;
			}
		}
	}
	return new MeshLib::Mesh("SubsurfaceMesh", new_nodes, new_elems);
}

int MshLayerMapper::LayerMapping(MeshLib::Mesh* new_mesh, const std::string &rasterfile,
                                 const unsigned nLayers, const unsigned layer_id, double noDataReplacementValue = 0.0)
{
	if (new_mesh == nullptr)
	{
		ERR("MshLayerMapper::LayerMapping() - Passed Mesh is NULL.");
		return 0;
	}

	if (nLayers >= layer_id)
	{
		GeoLib::Raster *raster(GeoLib::Raster::getRasterFromASCFile(rasterfile));
		if (! raster) {
			ERR("MshLayerMapper::LayerMapping - could not read raster file %s", rasterfile.c_str());
			return 0;
		}
		const double x0(raster->getOrigin()[0]);
		const double y0(raster->getOrigin()[1]);
		const double delta(raster->getRasterPixelDistance());
		const double no_data(raster->getNoDataValue());
		const std::size_t width(raster->getNCols());
		const std::size_t height(raster->getNRows());
		double const*const elevation(raster->begin());

		const std::pair<double, double> xDim(x0, x0 + width * delta); // extension in x-dimension
		const std::pair<double, double> yDim(y0, y0 + height * delta); // extension in y-dimension

		const size_t nNodes = new_mesh->getNNodes();
		const size_t nNodesPerLayer = nNodes / (nLayers+1);

		const size_t firstNode = layer_id * nNodesPerLayer;
		const size_t lastNode  = firstNode + nNodesPerLayer;

		const double half_delta = 0.5*delta;
		const std::vector<MeshLib::Node*> nodes = new_mesh->getNodes();
		for (unsigned i = firstNode; i < lastNode; ++i)
		{
			const double* coords = nodes[i]->getCoords();

			if (!isNodeOnRaster(*nodes[i], xDim, yDim))
			{
				if (layer_id == 0) // use default value
					nodes[i]->updateCoordinates(coords[0], coords[1], noDataReplacementValue);
				else // use z-value from layer above
					nodes[i]->updateCoordinates(coords[0], coords[1], (*nodes[i-nNodesPerLayer])[2]);
				continue;
			}

			// position in raster
			const double xPos ((coords[0] - xDim.first) / delta);
			const double yPos ((coords[1] - yDim.first) / delta);
			// raster cell index
			const size_t xIdx (static_cast<size_t>(floor(xPos)));
			const size_t yIdx (static_cast<size_t>(floor(yPos)));

			// deviation of mesh node from centre of raster cell ( in [-1:1) because it is normalised by delta/2 )
			const double xShift = (xPos-xIdx-half_delta)/half_delta;
			const double yShift = (yPos-yIdx-half_delta)/half_delta;

			const int xShiftIdx = static_cast<int>((xShift>=0) ? ceil(xShift) : floor(xShift));
			const int yShiftIdx = static_cast<int>((yShift>=0) ? ceil(yShift) : floor(yShift));

			// determining the neighbouring pixels that add weight to the interpolation
			const int x_nb[4] = {0, xShiftIdx, xShiftIdx, 0};
			const int y_nb[4] = {0, 0, yShiftIdx, yShiftIdx};

			double locZ[4];
			locZ[0] = elevation[yIdx*width + xIdx];
			if (fabs(locZ[0] - no_data) > std::numeric_limits<double>::min())
			{
				for (unsigned j=1; j<4; ++j)
				{
					locZ[j] = elevation[(yIdx+y_nb[j])*width + (xIdx+x_nb[j])];
					if (fabs(locZ[j] - no_data) < std::numeric_limits<double>::min())
						locZ[j]=locZ[0];
				}

				double ome[4];
				double xi = 1-fabs(xShift);
				double eta = 1-fabs(xShift);
				MathLib::MPhi2D(ome, xi, eta);

				double z(0.0);
				for(unsigned j = 0; j < 4; ++j)
					z += ome[j] * locZ[j];

				nodes[i]->updateCoordinates(coords[0], coords[1], z);
			}
			else
			{
				if (layer_id == 0) // use default value
					nodes[i]->updateCoordinates(coords[0], coords[1], noDataReplacementValue);
				else // use z-value from layer above
					nodes[i]->updateCoordinates(coords[0], coords[1], (*nodes[i-nNodesPerLayer])[2]);
			}
		}

		delete raster;
		return 1;
	}
	else
		ERR("MshLayerMapper::LayerMapping() - Mesh has only %d Layers, cannot assign layer %d.", nLayers, layer_id);
	return 0;
}

bool MshLayerMapper::isNodeOnRaster(const MeshLib::Node &node,
                                    const std::pair<double, double> &xDim,
                                    const std::pair<double, double> &yDim)
{
	if (node[0] < xDim.first || node[0] > xDim.second || node[1] < yDim.first || node[1] > yDim.second)
		return false;

	return true;
}

MeshLib::Mesh* MshLayerMapper::blendLayersWithSurface(MeshLib::Mesh* mesh, const unsigned nLayers, const std::string &dem_raster)
{
	// construct surface mesh from DEM
	const double dir[3] = {0,0,1};
	MeshLib::Mesh* dem = MeshLib::MeshSurfaceExtraction::getMeshSurface(*mesh, dir);
	MshLayerMapper::LayerMapping(dem, dem_raster, 0, 0);
	const std::vector<MeshLib::Node*> dem_nodes (dem->getNodes());

	const std::vector<MeshLib::Node*> mdl_nodes (mesh->getNodes());
	const size_t nNodes = mesh->getNNodes();
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
	//		else remove node (i.e. don't copy them)
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
	const std::vector<MeshLib::Element*> mdl_elements (mesh->getElements());
	const size_t nElems = mesh->getNElements();
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

			// construct pyrmid element based on missing node
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



