/**
 * \file MshLayerMapper.cpp
 * 01/11/2010 KR Initial implementation
 */

#include "MshLayerMapper.h"
#include "VtkRaster.h"

#include "MshEditor.h"
#include "matrix_class.h"
#include "msh_mesh.h"

#include <QImage>

MeshLib::CFEMesh* MshLayerMapper::CreateLayers(const MeshLib::CFEMesh* mesh,
                                               size_t nLayers,
                                               double thickness)
{
	if (nLayers < 1 || thickness <= 0)
	{
		std::cout <<
		"Error in MshLayerMapper::CreateLayers() - Invalid parameter: nLayers > 0 and thickness > 0 are required."
		          << std::endl;
		return NULL;
	}
/*
    if ((mesh->ele_vector[0]->GetElementType() != MshElemType::TRIANGLE) && (mesh->ele_vector[0]->GetElementType() != MshElemType::QUAD)) // check if mesh elements are triangles or quads
    {
        std::cout << "Error in MshLayerMapper::CreateLayers() - Method can only handle triangle- or quad-meshes... " << std::endl;
        return NULL;
    }
 */
	MeshLib::CFEMesh* new_mesh ( new MeshLib::CFEMesh() );
	const size_t nNodes = mesh->nod_vector.size();
	const size_t nElems = mesh->ele_vector.size();

	for (size_t layer_id = 0; layer_id < nLayers; layer_id++)
	{
		// add nodes for new layer
		size_t node_offset ( nNodes * layer_id );
		double z_offset ( layer_id*thickness );
		for (size_t i = 0; i < nNodes; i++)
		{
			const double* coords = mesh->nod_vector[i]->getData();
			new_mesh->nod_vector.push_back( new MeshLib::CNode(node_offset + i,
															   coords[0], 
															   coords[1], 
															   coords[2]-z_offset) );
		}

		if (layer_id > 0) // starting with the 2nd layer prism (or hex) elements can be created
		{
			// create prism elements connecting the last layer with the current one
			node_offset = (layer_id - 1) * nNodes;
			for (size_t i = 0; i < nElems; i++)
			{
				const MeshLib::CElem* sfc_elem( mesh->ele_vector[i] );
				MeshLib::CElem* elem( new MeshLib::CElem() );
				size_t nElemNodes = sfc_elem->getNodeIndices().Size();
				if (sfc_elem->GetElementType() == MshElemType::TRIANGLE)
					elem->setElementProperties(MshElemType::PRISM);                                           // extrude triangles to prism
				else if (sfc_elem->GetElementType() == MshElemType::QUAD)
					elem->setElementProperties(MshElemType::HEXAHEDRON);                                            // extrude quads to hexes
				else if (sfc_elem->GetElementType() == MshElemType::LINE)
					continue;                                            // line elements are ignored and not duplicated
				else
				{
					std::cout << "Error in MshLayerMapper::CreateLayers() - Method can only handle 2D mesh elements ..." << std::endl;
					std::cout << "Element " << i << " is of type \"" << MshElemType2String(sfc_elem->GetElementType()) << "\"." << std::endl;
					delete new_mesh;
					return NULL;
				}
				elem->SetPatchIndex(layer_id - 1);
				elem->SetNodesNumber(2 * nElemNodes);
				elem->getNodeIndices().resize(2 * nElemNodes);
				for (size_t j = 0; j < nElemNodes; j++)
				{
					long idx = sfc_elem->GetNodeIndex(j);
					elem->SetNodeIndex(j, node_offset + idx);
					elem->SetNodeIndex(j + nElemNodes, node_offset + nNodes + idx);
				}
				new_mesh->ele_vector.push_back(elem);
			}
		}
	}

	new_mesh->setNumberOfNodesFromNodesVectorSize ();
	new_mesh->setNumberOfMeshLayers(nLayers);

	new_mesh->ConstructGrid();
	new_mesh->FillTransformMatrix();

	return new_mesh;
}

MeshLib::CFEMesh* MshLayerMapper::LayerMapping(const MeshLib::CFEMesh* msh,
                                               const std::string &rasterfile,
                                               const size_t nLayers,
                                               const size_t layer_id,
                                               bool removeNoDataValues)
{
	if (msh == NULL)
		return NULL;
	if (msh->getNumberOfMeshLayers() >= layer_id)
	{
		if (msh == NULL)
		{
			std::cout <<
			"Error in MshLayerMapper::LayerMapping() - Passed Mesh is NULL..." <<
			std::endl;
			return NULL;
		}
		MeshLib::CFEMesh* new_mesh( new MeshLib::CFEMesh(*msh) );

		double x0(0), y0(0), delta(1);
		size_t width(1), height(1);
		float* elevation = VtkRaster::loadDataFromASC(rasterfile, x0, y0, width,height, delta);

		if (elevation == NULL)
		{
			delete [] elevation;
			return NULL;
		}

		std::pair<double, double> xDim(x0, x0 + width * delta); // extension in x-dimension
		std::pair<double, double> yDim(y0, y0 + height * delta); // extension in y-dimension

		if (!meshFitsImage(new_mesh, xDim, yDim))
		{
			delete [] elevation;
			return NULL;
		}

		size_t nNodes = msh->nod_vector.size();
		size_t nNodesPerLayer = nNodes / nLayers;

		size_t firstNode = layer_id * nNodesPerLayer;
		size_t lastNode  = firstNode + nNodesPerLayer;

		std::vector<size_t> noData_nodes;
		const double half_delta = 0.5*delta;
		for(size_t i = firstNode; i < lastNode; i++)
		{
			const double* coords = msh->nod_vector[i]->getData();
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
			locZ[0] = elevation[2*(yIdx*width + xIdx)];
			if (fabs(locZ[0] + 9999) > std::numeric_limits<double>::min())
			{
				for (size_t j=1; j<4; j++)
				{
					locZ[j] = elevation[2*((yIdx+y_nb[j])*width + (xIdx+x_nb[j]))];
					if (fabs(locZ[j] + 9999) < std::numeric_limits<double>::min())
						locZ[j]=locZ[0];
				}

				double ome[4];
				double xi = 1-fabs(xShift);
				double eta = 1-fabs(xShift);
				MPhi2D(ome, xi, eta);

				double z(0.0);
				for(size_t j = 0; j < 4; j++)
					z += ome[j] * locZ[j];
				new_mesh->nod_vector[i]->SetZ(z);
				new_mesh->nod_vector[i]->SetMark(true);
			}
			else
			{
				new_mesh->nod_vector[i]->SetZ(0);
				new_mesh->nod_vector[i]->SetMark(false);
				noData_nodes.push_back(i);
			}
		}

		if ((nLayers == 1) && removeNoDataValues)
		{
			if (noData_nodes.size() < (new_mesh->nod_vector.size() - 2))
			{
				std::cout << "Warning: Removing " << noData_nodes.size() <<
				" mesh nodes at NoData values ... " << std::endl;
				MeshLib::CFEMesh* red_mesh = MshEditor::removeMeshNodes(
				        new_mesh,
				        noData_nodes);
				if (!new_mesh->ele_vector.empty())
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

		new_mesh->ConstructGrid();
		new_mesh->FillTransformMatrix();

		delete [] elevation;
		return new_mesh;
	}
	else
		std::cout << "Error in MshLayerMapper::LayerMapping() - Mesh has only " <<
		msh->getNumberOfMeshLayers() << " Layers, cannot assign layer " << layer_id <<
		"..." << std::endl;
	return NULL;
}

// KR, based on code by WW (Note: this method has not been tested yet and will probably fail miserably!)
bool MshLayerMapper::meshFitsImage(const MeshLib::CFEMesh* msh,
                                   const std::pair<double, double> &xDim,
                                   const std::pair<double, double> &yDim)
{
	double const* pnt (msh->nod_vector[0]->getData());
	double xMin(pnt[0]);
	double yMin(pnt[1]);
	double xMax(pnt[0]);
	double yMax(pnt[1]);

	size_t nNodes = msh->nod_vector.size();
	for (size_t i = 1; i < nNodes; i++)
	{
		pnt = msh->nod_vector[i]->getData();
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
		std::cout << "Warning: Mesh does not fit into given raster file." << std::endl;
		return false;
	}
	return true;
}

void MshLayerMapper::CheckLayerMapping(MeshLib::CFEMesh* mesh, const size_t nLayers, int integ)
{
	double ref_dep = -999999999.0;

	size_t nNodesPerLayer = mesh->nod_vector.size() / (nLayers + 1);

	//18.02.2009 WW
	if (integ)
	{
		for (size_t i = 0; i < nNodesPerLayer; i++)
			for (size_t k = 0; k < nLayers; k++)
			{
				MeshLib::CNode* node = mesh->nod_vector[k * nNodesPerLayer + i];
				if (k == 0)
					node->SetBoundaryType('0');  // on the surface
				else if (k == (nLayers - 1))
					node->SetBoundaryType('1');  // on the bottom
				else
					node->SetBoundaryType('I');  // interior node
			}
	}

	size_t flat(0);
	for (size_t i = 0; i < nNodesPerLayer; i++)
	{
		std::vector<long> tmp_connected_nodes;
		flat = 0;

		for (size_t k = 0; k < nLayers - 2; k++) // top layer is not checked
		{
			MeshLib::CNode* bNode = mesh->nod_vector[k * nNodesPerLayer + i]; // current node
			MeshLib::CNode* tNode = mesh->nod_vector[(k + 1) * nNodesPerLayer + i]; // same node but one layer below

			if (!tNode->GetMark())
			{
				if (k == 0)
				{
					tmp_connected_nodes.clear();
					std::vector<size_t> const& connected_nodes (
					        tNode->getConnectedNodes());
					const size_t n_connected_nodes (connected_nodes.size());
					for (size_t j = 0; j < n_connected_nodes; j++)
						tmp_connected_nodes.push_back(connected_nodes[j]);
				}

				tNode->SetZ(bNode->getData()[2]); // z coordinate changed to layer above
				tNode->getConnectedNodes().clear();
				for (int j = k; j >= 0; j--) //WW/YW. 23.01.2009
				{
					MeshLib::CNode* nNode =
					        mesh->nod_vector[j * nNodesPerLayer + i];
					if (nNode->GetMark())
					{
						tNode->getConnectedNodes().push_back(
						        nNode->GetIndex());
						break;
					}
				}
				flat++;
			}
		}

		//---- 27.01.2009. WW
		if (flat == nLayers - 2 /*1*/)
		{
			MeshLib::CNode* bNode = mesh->nod_vector[nNodesPerLayer + i];
			bNode->SetMark(true);
			bNode->getConnectedNodes().clear();
			for (size_t j = 0; j < tmp_connected_nodes.size(); j++)
				bNode->getConnectedNodes().push_back(tmp_connected_nodes[j]);

			MeshLib::CNode* tNode =
			        mesh->nod_vector[(nLayers - 1) * nNodesPerLayer + i];
			tNode->SetMark(false);
			bNode->SetZ(tNode->getData()[2]);
			bNode->SetBoundaryType('1');
			//
			for (size_t k = 1; k < nLayers; k++)
			{
				tNode = mesh->nod_vector[(k + 1) * nNodesPerLayer + i];
				tNode->SetZ(ref_dep);
				tNode->getConnectedNodes().clear();
				tNode->getConnectedNodes().push_back(bNode->GetIndex());
			}
		}
	}

	std::vector<MeshLib::CElem*> new_elems;
	std::vector<size_t> false_node_idx;
	size_t nElems = mesh->ele_vector.size();
	for(size_t i = 0; i < nElems; i++)
	{
		MeshLib::CElem* elem = mesh->ele_vector[i];
		elem->SetMark(true);

		flat = 0;
		for (size_t k = 0; k < 3; k++) //skip element when one node is marked ref_dep

			if ( fabs(elem->GetNode(k)->getData()[2]  + ref_dep) <
			     std::numeric_limits<double>::min() ||
			     fabs(elem->GetNode(k + 3)->getData()[2] + ref_dep) <
			     std::numeric_limits<double>::min() )
			{
				flat = 1;
				elem->SetMark(false);
				break;
			}
		if (flat == 1)
			continue;

		// If all nodes are okay, check if two z-values are identical
		for (size_t k = 0; k < 3; k++)
			if(fabs(elem->GetNode(k +
			                      3)->getData()[2] - elem->GetNode(k)->getData()[2]) <
			   std::numeric_limits<double>::min())
				false_node_idx.push_back(k);

		switch(false_node_idx.size())
		{
		case 0: // everything okay, take the prism as it is
		{
			elem->SetMark(true);
			break;
		}
		case 1: // one node of the prism is marked false, i.e. create two tetrahedron elements from the remaining 5 prism nodes
		{
			size_t a = false_node_idx[0];
			size_t b = (a + 2) % 3;
			size_t c = (a + 1) % 3;

			if (elem->GetNode(a + 3)->GetBoundaryType() == '1') //24.02.2009. WW
				elem->GetNode(a)->SetBoundaryType('1');

			// create a new tetrahedron
			MeshLib::CElem* new_elem( new MeshLib::CElem() );
			new_elem->SetMark(true);
			new_elem->SetElementType(MshElemType::TETRAHEDRON);
			new_elem->SetPatchIndex(elem->GetPatchIndex());
			new_elem->SetBoundaryType('I');
			new_elem->SetNodesNumber(4);

			Math_Group::vec<MeshLib::CNode*> nodes(4);
			nodes[0] = mesh->nod_vector[elem->getNodeIndices()[a]];
			nodes[1] = mesh->nod_vector[elem->getNodeIndices()[b + 3]];
			nodes[2] = mesh->nod_vector[elem->getNodeIndices()[c + 3]];
			nodes[3] = mesh->nod_vector[elem->getNodeIndices()[c]];
			new_elem->SetNodes(nodes, true);

			new_elem->getNodeIndices().resize(4);
			for (size_t k = 0; k < 4; k++)
				new_elem->getNodeIndices()[k] = elem->GetNode(k)->GetIndex();
			new_elems.push_back(new_elem);

			// change prism-element to 2nd tetrahedron
			elem->SetMark(true);
			elem->SetElementType(MshElemType::TETRAHEDRON);
			elem->SetNodesNumber(4);

			nodes[0] = mesh->nod_vector[elem->getNodeIndices()[b]];
			nodes[1] = mesh->nod_vector[elem->getNodeIndices()[a]];
			nodes[2] = mesh->nod_vector[elem->getNodeIndices()[c]];
			nodes[3] = mesh->nod_vector[elem->getNodeIndices()[b + 3]];
			elem->SetNodes(nodes, true);

			elem->getNodeIndices().resize(4);
			for (size_t k = 0; k < 4; k++)
				elem->getNodeIndices()[k] = elem->GetNode(k)->GetIndex();
			break;
		}
		case 2: // two nodes of the prism are marked false, i.e. create a tetrahedron element from the remaining 4 prism nodes
		{
			size_t a = false_node_idx[0];
			size_t b = (a + 2) % 3;
			size_t c = (a + 1) % 3;
			if (false_node_idx[1] == b)
				a = c;
			else if(false_node_idx[1] == c)
				a = b;

			elem->SetMark(true);
			elem->SetElementType(MshElemType::TETRAHEDRON);
			elem->SetNodesNumber(4);

			Math_Group::vec<MeshLib::CNode*> nodes(4);
			nodes[0] = mesh->nod_vector[elem->getNodeIndices()[b]];
			nodes[1] = mesh->nod_vector[elem->getNodeIndices()[a]];
			nodes[2] = mesh->nod_vector[elem->getNodeIndices()[c]];
			nodes[3] = mesh->nod_vector[elem->getNodeIndices()[a + 3]];
			elem->SetNodes(nodes, true);

			elem->getNodeIndices().resize(4);
			for (size_t k = 0; k < 4; k++)
				elem->getNodeIndices()[k] = elem->GetNode(k)->GetIndex();
			/*
			   //for j, l nodes if they becomes on top surface. 24.02.2009. WW
			   if (node_b->GetBoundaryType()=='1')
			    elem->nodes[0]->SetBoundaryType('1');
			   if (node_t->GetBoundaryType()=='1')
			    elem->nodes[2]->SetBoundaryType('1');
			 */
			break;
		}
		case 3: // three nodes of the prism is marked false, ditch the whole element
		{
			elem->SetMark(false);
			break;
		}
		}
	}

	// add the newly created elements to the elements vector
	for(size_t i = 0; i < new_elems.size(); i++)
		mesh->ele_vector.push_back(new_elems[i]);

	// correct indeces of elements and delete false elements
	std::vector<MeshLib::CElem*>::iterator beg_e = mesh->ele_vector.begin( ), last_e;
	long counter = 0;
	while ( beg_e != mesh->ele_vector.end() )
	{
		last_e = beg_e++;
		MeshLib::CElem* elem = *last_e;
		if (elem->GetMark())
		{
			elem->SetIndex(counter);
			counter++;
			/* KR unused variable
			   for (int j=0; j<elem->GetVertexNumber(); j++)
			   {
			    if (!elem->GetNode(j)->GetMark())
			    {
			        MeshLib::CNode* node = mesh->nod_vector[elem->GetNode(j)->connected_nodes[0]];
			    }
			   }
			 */
		}
		else
		{
			delete elem;
			beg_e = mesh->ele_vector.erase(last_e);
		}
	}

	// correct indeces of nodes and delete false nodes
	counter = 0;
	std::vector<MeshLib::CNode*>::iterator beg = mesh->nod_vector.begin( ), last;
	while ( beg != mesh->nod_vector.end( ) )
	{
		last = beg++;
		MeshLib::CNode* node = *last;
		if (node->GetMark())
		{
			node->SetIndex(counter);
			node->getConnectedElementIDs().clear();
			node->getConnectedNodes().clear();
			counter++;
		}
		else
		{
			delete  node;
			node = NULL;
			beg = mesh->nod_vector.erase(last);
		}
	}

	// correct element indeces again after deleting nodes
	nElems = mesh->ele_vector.size();
	for(size_t i = 0; i < nElems; i++)
		for(int j = 0; j < mesh->ele_vector[i]->GetVertexNumber(); j++)
			mesh->ele_vector[i]->getNodeIndices()[j] =
			        mesh->ele_vector[i]->GetNode(j)->GetIndex();

	// delete elements if two of its nodes are identical (can this actually happen?????)
	beg_e = mesh->ele_vector.begin();
	counter = 0;
	bool flatf = false;
	while ( beg_e != mesh->ele_vector.end( ) )
	{
		last_e = beg_e++;
		MeshLib::CElem* elem = *last_e;

		//10.02.2009. WW !!!!!!!!!!!!!!!!!!!!!!
		for (int j = 0; j < elem->GetVertexNumber(); j++)
		{
			flatf = false;
			for (int k = j; k < elem->GetVertexNumber(); k++)
				if (elem->GetNodeIndex(j) == elem->GetNodeIndex(k))
				{
					flatf = true;
					break;
				}
		}
		if (flatf)
		{
			delete elem;
			beg_e = mesh->ele_vector.erase(last_e);
		}
		else
		{
			elem->SetIndex(counter);
			counter++;
			/* KR unused variable
			   for (int j=0; j<elem->GetVertexNumber(); j++)
			   {
			   if (!elem->GetNode(j)->GetMark())
			   {
			       MeshLib::CNode* node = mesh->nod_vector[elem->GetNode(j)->connected_nodes[0]];
			   }
			   }
			 */
		}
	}

	mesh->ConnectedElements2Node();

	// delete nodes that are not connected to any element (can this happen???)
	counter = 0;
	beg = mesh->nod_vector.begin( );
	while ( beg != mesh->nod_vector.end( ) )
	{
		last = beg++;
		MeshLib::CNode* node = *last;
		if ( node->getConnectedElementIDs().empty() )
		{
			delete node;
			node = NULL;
			beg = mesh->nod_vector.erase(last);
		}
		else
		{
			node->SetIndex(counter);
			counter++;
		}
	}

	nElems = mesh->ele_vector.size();
	for (size_t i = 0; i < nElems; i++)
		for (int j = 0; j < mesh->ele_vector[i]->GetVertexNumber(); j++)
			mesh->ele_vector[i]->getNodeIndices()[j] =
			        mesh->ele_vector[i]->GetNode(j)->GetIndex();

	mesh->ConstructGrid();
}
