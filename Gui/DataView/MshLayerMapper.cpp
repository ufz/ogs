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
#include "MshEditor.h"
#include "MathTools.h"

#include <QImage>

MeshLib::Mesh* MshLayerMapper::CreateLayers(const MeshLib::Mesh* mesh, unsigned nLayers, float* thickness)
{
	if (nLayers < 1 || thickness <= 0 || mesh->getDimension() != 2)
	{
		std::cout << "Error in MshLayerMapper::CreateLayers() - A 2D mesh with nLayers > 0 is required as input." << std::endl;
		return NULL;
	}

	const size_t nNodes = mesh->getNNodes();
	const size_t nElems = mesh->getNElements();
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

				for (unsigned i = 0; i < nElems; ++i)
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
						if (sfc_elem->getGeomType() == MshElemType::TRIANGLE)	// extrude triangles to prism
							new_elems[elem_offset+i] = new MeshLib::Prism(e_nodes, mat_id);
						else if (sfc_elem->getGeomType() == MshElemType::QUAD)	// extrude quads to hexes
							new_elems[elem_offset+i] = new MeshLib::Hex(e_nodes, mat_id);
					}
					else
					{
						std::cout << "Warning in MshLayerMapper::CreateLayers() - Method can only handle 2D mesh elements ..." << std::endl;
						std::cout << "Skipping Element " << i << " of type \"" << MshElemType2String(sfc_elem->getGeomType()) << "\"." << std::endl;
					}
				}
			}
			else
			{
				std::cout << "Error in MshLayerMapper::CreateLayers() - Layer thickness for layer "
					      << (layer_id-1) << " is " << thickness[layer_id-1] << " (needs to be >0)." << std::endl;
				return NULL;
			}
		}
	}
	return new MeshLib::Mesh("SubsurfaceMesh", new_nodes, new_elems);
}

int MshLayerMapper::LayerMapping(MeshLib::Mesh* new_mesh, const std::string &rasterfile,
                                 const unsigned nLayers, const unsigned layer_id, bool removeNoDataValues)
{
	if (new_mesh == NULL)
	{
		std::cout <<
		"Error in MshLayerMapper::LayerMapping() - Passed Mesh is NULL..." <<
		std::endl;
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

		if (!meshFitsImage(new_mesh, xDim, yDim))
		{
			delete [] elevation;
			return 0;
		}

		const size_t nNodes = new_mesh->getNNodes();
		const size_t nNodesPerLayer = nNodes / (nLayers+1);

		const size_t firstNode = layer_id * nNodesPerLayer;
		const size_t lastNode  = firstNode + nNodesPerLayer;

		std::vector<size_t> noData_nodes;
		const double half_delta = 0.5*delta;
		const std::vector<MeshLib::Node*> nodes = new_mesh->getNodes();
		for (unsigned i = firstNode; i < lastNode; ++i)
		{
			const double* coords = nodes[i]->getCoords();
			// position in raster
			const double xPos ((coords[0] - xDim.first) / delta);
			const double yPos ((coords[1] - yDim.first) / delta);
			// raster cell index
			const size_t xIdx (static_cast<size_t>(floor(xPos)));
			const size_t yIdx (static_cast<size_t>(floor(yPos)));

			// deviation of mesh node from centre of raster cell ( in [-1:1) because it is normalised by delta/2 )
			const double xShift = (xPos-xIdx-half_delta)/half_delta;
			const double yShift = (yPos-yIdx-half_delta)/half_delta;

			const int xShiftIdx = (xShift>=0) ? ceil(xShift) : floor(xShift);
			const int yShiftIdx = (yShift>=0) ? ceil(yShift) : floor(yShift);

			// determining the neighbouring pixels that add weight to the interpolation
			const size_t x_nb[4] = {0, xShiftIdx, xShiftIdx, 0};
			const size_t y_nb[4] = {0, 0, yShiftIdx, yShiftIdx};

			double locZ[4];
			locZ[0] = elevation[yIdx*width + xIdx];
			if (fabs(locZ[0] + no_data) > std::numeric_limits<double>::min())
			{
				for (unsigned j=1; j<4; ++j)
				{
					locZ[j] = elevation[(yIdx+y_nb[j])*width + (xIdx+x_nb[j])];
					if (fabs(locZ[j] + no_data) < std::numeric_limits<double>::min())
						locZ[j]=locZ[0];
				}

				double ome[4];
				double xi = 1-fabs(xShift);
				double eta = 1-fabs(xShift);
				MathLib::MPhi2D(ome, xi, eta);

				double z(0.0);
				for(unsigned j = 0; j < 4; ++j)
					z += ome[j] * locZ[j];
				const double* coords (nodes[i]->getCoords());
				nodes[i]->updateCoordinates(coords[0], coords[1], z);
			}
			else
			{
				const double* coords (nodes[i]->getCoords());
				nodes[i]->updateCoordinates(coords[0], coords[1], 0);
				noData_nodes.push_back(i);
			}
		}

		if ((nLayers == 0) && removeNoDataValues)
		{
			if (noData_nodes.size() < (nNodes - 2))
			{
				std::cout << "Warning: Removing " << noData_nodes.size()
					      << " mesh nodes at NoData values ... " << std::endl;
				MeshLib::Mesh* red_mesh = MeshLib::MshEditor::removeMeshNodes(new_mesh, noData_nodes);
				if (new_mesh->getNElements() == 0)
				{
					delete new_mesh;
					new_mesh = red_mesh;
				}
				else
				{
					delete red_mesh;
					std::cout << "Too many NoData values..." << std::endl;
				}
			}
			else
				std::cout << "Too many NoData values..." << std::endl;
		}

		delete raster;
		return 1;
	}
	else
		std::cout << "Error in MshLayerMapper::LayerMapping() - Mesh has only "
		          << nLayers << " Layers, cannot assign layer " << layer_id
				  << "..." << std::endl;
	return 0;
}

bool MshLayerMapper::meshFitsImage(const MeshLib::Mesh* msh,
                                   const std::pair<double, double> &xDim,
                                   const std::pair<double, double> &yDim)
{
	const size_t nNodes = msh->getNNodes();
	const std::vector<MeshLib::Node*> nodes = msh->getNodes();
	const double* pnt;
	double xMin(std::numeric_limits<double>::max());
	double yMin(std::numeric_limits<double>::max());
	double xMax(std::numeric_limits<double>::min());
	double yMax(std::numeric_limits<double>::min());

	for (unsigned i = 1; i < nNodes; ++i)
	{
		pnt = nodes[i]->getCoords();
		if (xMin > pnt[0])
			xMin = pnt[0];
		else if (xMax < pnt[0])
			xMax = pnt[0];

		if (yMin > pnt[1])
			yMin = pnt[1];
		else if (yMax < pnt[1])
			yMax = pnt[1];
	}

	if (xMin < xDim.first || xMax > xDim.second || yMin < yDim.first || yMax > yDim.second)
	{
		std::cout << "Warning: Extension of mesh is larger than extension of given raster file." << std::endl;
		std::cout << "Mesh Extend: (" << xMin << ", " << yMin << "):(" << xMax << ", " << yMax << ")" << std::endl;
		std::cout << "Raster Extend: (" << xDim.first << ", " << yDim.first << "):(" << xDim.second << ", " << yDim.second << ")" << std::endl;
		return false;
	}
	return true;
}

MeshLib::Mesh* MshLayerMapper::blendLayersWithSurface(MeshLib::Mesh* mesh, const unsigned nLayers, const std::string &dem_raster)
{
	// construct surface mesh from DEM
	const double dir[3] = {0,0,1};
	MeshLib::Mesh* dem = MeshLib::MshEditor::getMeshSurface(*mesh, dir);
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
			std::cout << "Warning: Node " << i << " (in bottom-layer) is above surface node " << (i-bottom_firstNode) << ". (" << coords[2] << " > " << dem_coords[2] << ")" << std::endl;
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



