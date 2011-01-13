/**
 * \file MshLayerMapper.cpp
 * 01/11/2010 KR Initial implementation
 */

#include "MshLayerMapper.h"
#include "OGSRaster.h"

#include "msh_mesh.h"
#include "matrix_class.h"

#include <QImage>


Mesh_Group::CFEMesh* MshLayerMapper::CreateLayers(const Mesh_Group::CFEMesh* mesh, size_t nLayers, double thickness)
{
	if (nLayers < 1 || thickness <= 0)
	{
		std::cout << "Error in MshLayerMapper::CreateLayers() - Invalid parameter: nLayers > 0 and thickness > 0 are required." << std::endl;
		return NULL;
	}

	if ((mesh->ele_vector[0]->GetElementType() != MshElemType::TRIANGLE) && (mesh->ele_vector[0]->GetElementType() != MshElemType::QUAD)) // check if mesh elements are triangles or quads
	{
		std::cout << "Error in MshLayerMapper::CreateLayers() - Method can only handle triangle- or quad-meshes... " << std::endl;
		return NULL;
	}

	Mesh_Group::CFEMesh* new_mesh ( new Mesh_Group::CFEMesh() );
	size_t nNodes = mesh->nod_vector.size();
	size_t nElems = mesh->ele_vector.size();
	size_t nElemNodes = mesh->ele_vector[0]->nodes_index.Size();

	for (size_t layer_id=0; layer_id<nLayers; layer_id++)
	{
		// add nodes for new layer
		size_t node_offset ( nNodes*layer_id );
		double z_offset ( layer_id*thickness );
		for (size_t i=0; i<nNodes; i++)
		{
			Mesh_Group::CNode* node = new Mesh_Group::CNode( node_offset + i );
			double coords[3] = { mesh->nod_vector[i]->X(), mesh->nod_vector[i]->Y(), mesh->nod_vector[i]->Z()-z_offset };
			node->SetCoordinates(coords);
			new_mesh->nod_vector.push_back(node);
		}

		if (layer_id>0) // starting with the 2nd layer prism (or hex) elements can be created
		{
			// create prism elements connecting the last layer with the current one
			node_offset = (layer_id-1)*nNodes;
			for (size_t i=0; i<nElems; i++)
			{
				Mesh_Group::CElem* elem = new Mesh_Group::CElem();
				if (mesh->ele_vector[i]->GetElementType()==MshElemType::TRIANGLE) elem->SetElementType(MshElemType::PRISM); // extrude triangles to prism
				else elem->SetElementType(MshElemType::HEXAHEDRON); // extrude quads to hexes
				elem->SetPatchIndex(layer_id-1);
				elem->SetNodesNumber(2*nElemNodes);
				elem->nodes_index.resize(2*nElemNodes);
				for (size_t j=0; j<nElemNodes; j++)
				{
					long idx = mesh->ele_vector[i]->GetNodeIndex(j);
					elem->SetNodeIndex(j, node_offset+idx);
					elem->SetNodeIndex(j+nElemNodes, node_offset+nNodes+idx);
				}
				new_mesh->ele_vector.push_back(elem);
			}
		}
	}

	new_mesh->setNumberOfMeshLayers(nLayers);

	new_mesh->ConstructGrid();
	new_mesh->FillTransformMatrix();

	return new_mesh;
}

// KR, based on code by WW
Mesh_Group::CFEMesh* MshLayerMapper::LayerMapping(const Mesh_Group::CFEMesh* msh, const std::string &rasterfile, const size_t nLayers, const size_t layer_id)
{
	if (msh == NULL) return NULL;
	if (msh->getNumberOfMeshLayers() >= layer_id)
	{
		if (msh==NULL)
		{
			std::cout << "Error in MshLayerMapper::LayerMapping() - Passed Mesh is NULL..." << std::endl;
			return NULL;
		}
		Mesh_Group::CFEMesh* new_mesh = new Mesh_Group::CFEMesh(*msh);

		double x0(0), y0(0), delta(1);
		size_t width(1), height(1);
		double* elevation = OGSRaster::loadDataFromASC(QString::fromStdString(rasterfile), x0, y0, width, height, delta);

		if (elevation==NULL)
		{
			delete []elevation;
			return NULL;
		}

		std::pair<double, double> xDim(x0, x0+width*delta);		// extension in x-dimension
		std::pair<double, double> yDim(y0, y0+height*delta);	// extension in y-dimension

		if (!meshFitsImage(msh, xDim, yDim))
		{
			delete []elevation;
			return NULL;
		}

		double locX[4];
		double locY[4];
		double locZ[4];

		size_t nNodes = msh->nod_vector.size();
		size_t nNodesPerLayer = nNodes / nLayers;

		size_t firstNode = layer_id * nNodesPerLayer;
		size_t lastNode  = firstNode + nNodesPerLayer;

		for(size_t i=firstNode; i<lastNode; i++)
		{
			size_t xPos (static_cast<size_t>(ceil((msh->nod_vector[i]->X() - xDim.first) / delta)));
			size_t yPos (static_cast<size_t>(ceil((msh->nod_vector[i]->Y() - yDim.first) / delta)));

			locX[0] = xDim.first+xPos*delta;
			locY[0] = yDim.first+yPos*delta;
			locZ[0] = elevation[yPos*width+xPos];

			locX[1] = xDim.first+(xPos+1)*delta;
			locY[1] = yDim.first+yPos*delta;
			locZ[1] = elevation[yPos*width+(xPos+1)];

			locX[2] = xDim.first+(xPos+1)*delta;
			locY[2] = yDim.first+(yPos+1)*delta;
			locZ[2] = elevation[(yPos+1)*width+(xPos+1)];

			locX[3] = xDim.first+xPos*delta;
			locY[3] = yDim.first+(yPos+1)*delta;
			locZ[3] = elevation[(yPos+1)*width+xPos];

			bool noData(false);
			for(size_t j=0; j<4; j++)
			{
				if(fabs(locZ[j]+9999)<std::numeric_limits<double>::min()) noData = true;
			}

			if(!noData)
			{
				 // Interpolate
				double ome[4];
				double xi  = 2.0*(msh->nod_vector[i]->X()-0.5*(locX[0]+locX[1]))/delta;
				double eta = 2.0*(msh->nod_vector[i]->Y()-0.5*(locY[1]+locY[2]))/delta;
				MPhi2D(ome, xi, eta);

				double z(0.0);
				for(size_t j=0; j<4; j++) z += ome[j]*locZ[j];
				new_mesh->nod_vector[i]->SetZ(z);
				new_mesh->nod_vector[i]->SetMark(true);
			}
			else
			{
				std::cout << "Warning: For node " << i << " (" << msh->nod_vector[i]->X() << ", " << msh->nod_vector[i]->Y() << ") there is no elevation data in the given raster." << std::endl;
				new_mesh->nod_vector[i]->SetZ(0);
				new_mesh->nod_vector[i]->SetMark(false);
			}
		}

		new_mesh->ConstructGrid();
		new_mesh->FillTransformMatrix();
		delete []elevation;
		return new_mesh;
	}
	else
	{
		std::cout << "Error in MshLayerMapper::LayerMapping() - Mesh has only " << msh->getNumberOfMeshLayers() << " Layers, cannot assign layer " << layer_id << "..." << std::endl;
	}
	return NULL;
}

// KR, based on code by WW (Note: this method has not been tested yet and will probably fail miserably!)
bool MshLayerMapper::meshFitsImage(const Mesh_Group::CFEMesh* msh, const std::pair<double, double> &xDim, const std::pair<double, double> &yDim)
{
	double xMin(msh->nod_vector[0]->X());
	double yMin(msh->nod_vector[0]->Y());
	double xMax(msh->nod_vector[0]->X());
	double yMax(msh->nod_vector[0]->Y());

	size_t nNodes = msh->nod_vector.size();
	for (size_t i=1; i<nNodes; i++)
	{
		if (xMin > msh->nod_vector[i]->X()) xMin = msh->nod_vector[i]->X();
		else if (xMax < msh->nod_vector[i]->X()) xMax = msh->nod_vector[i]->X();

		if (yMin > msh->nod_vector[i]->Y()) yMin = msh->nod_vector[i]->Y();
		else if (yMax < msh->nod_vector[i]->Y()) yMax = msh->nod_vector[i]->Y();
	}

	if (xMin < xDim.first || xMax > xDim.second || yMin < yDim.first || yMax > yDim.second)
	{
		std::cout << "Warning: Mesh does not fit into given raster file." << std::endl;
		return false;
	}
	return true;
}

void MshLayerMapper::CheckLayerMapping(Mesh_Group::CFEMesh* mesh, const size_t nLayers, int integ)
{
	double ref_dep = -999999999.0;

	size_t nNodesPerLayer = mesh->nod_vector.size() / (nLayers+1);

	//18.02.2009 WW
	if (integ)
	{
		for (size_t i=0; i<nNodesPerLayer; i++)
		{
		   for (size_t k=0; k<nLayers; k++)
		   {
			  Mesh_Group::CNode* node = mesh->nod_vector[k*nNodesPerLayer+i];
			  if (k==0) node->SetBoundaryType('0');	               // on the surface
			  else if (k==(nLayers-1)) node->SetBoundaryType('1'); // on the bottom
			  else node->SetBoundaryType('I');                     // interior node
		   }
		}
	}


	size_t flat(0);
	for (size_t i=0; i<nNodesPerLayer; i++)
	{
		std::vector<long> tmp_connected_nodes;
		flat = 0;

		for (size_t k=0; k<nLayers-2; k++) // top layer is not checked
		{
			Mesh_Group::CNode* bNode = mesh->nod_vector[k*nNodesPerLayer+i];     // current node
			Mesh_Group::CNode* tNode = mesh->nod_vector[(k+1)*nNodesPerLayer+i]; // same node but one layer below

			if (!tNode->GetMark())
			{
				if (k==0)
				{
					tmp_connected_nodes.clear();
					for (size_t j=0; j<tNode->connected_nodes.size(); j++)
						tmp_connected_nodes.push_back(tNode->connected_nodes[j]);
				}

				tNode->SetZ(bNode->Z()); // z coordinate changed to layer above
				tNode->connected_nodes.clear();
				for (int j=k; j>=0; j--)  //WW/YW. 23.01.2009
				{
					Mesh_Group::CNode* nNode = mesh->nod_vector[j*nNodesPerLayer+i];
					if (nNode->GetMark())
					{
						tNode->connected_nodes.push_back(nNode->GetIndex());
						break;
					}
				}
				flat++;
			}
		}

		//---- 27.01.2009. WW
		if (flat==nLayers-2/*1*/)
		{

			Mesh_Group::CNode* bNode = mesh->nod_vector[nNodesPerLayer+i];
			bNode->SetMark(true);
			bNode->connected_nodes.clear();
			for (size_t j=0; j<tmp_connected_nodes.size(); j++)
				bNode->connected_nodes.push_back(tmp_connected_nodes[j]);

			Mesh_Group::CNode* tNode = mesh->nod_vector[(nLayers-1)*nNodesPerLayer+i];
			tNode->SetMark(false);
			bNode->SetZ(tNode->Z());
			bNode->SetBoundaryType('1');
			//
			for (size_t k=1; k<nLayers; k++)
			{
				tNode = mesh->nod_vector[(k+1)*nNodesPerLayer+i];
				tNode->SetZ(ref_dep);
				tNode->connected_nodes.clear();
				tNode->connected_nodes.push_back(bNode->GetIndex());
			}
		}

	}


	std::vector<Mesh_Group::CElem*> new_elems;
	std::vector<size_t> false_node_idx;
	size_t nElems = mesh->ele_vector.size();
	for(size_t i=0; i<nElems; i++)
	{
		Mesh_Group::CElem* elem = mesh->ele_vector[i];
		elem->SetMark(true);

		flat = 0;
		for (size_t k=0; k<3;k++) //skip element when one node is marked ref_dep
		{
			if ( fabs(elem->GetNode(k)->Z()  +ref_dep)<std::numeric_limits<double>::min() ||
			     fabs(elem->GetNode(k+3)->Z()+ref_dep)<std::numeric_limits<double>::min() )
			{
				flat = 1;
				elem->SetMark(false);
				break;
			}
		}
		if (flat==1) continue;

		// If all nodes are okay, check if two z-values are identical
		for (size_t k=0; k<3;k++)
		{
			if(fabs(elem->GetNode(k+3)->Z()-elem->GetNode(k)->Z()) < std::numeric_limits<double>::min())
				false_node_idx.push_back(k);
		}

		switch(false_node_idx.size())
		{
			case 0: // everything okay, take the prism as it is
			{
				elem->SetMark(true);
				break;
			}
			case 1:	// one node of the prism is marked false, i.e. create two tetrahedron elements from the remaining 5 prism nodes
			{
				size_t a = false_node_idx[0];
				size_t b = (a+2)%3;
				size_t c = (a+1)%3;

				if (elem->GetNode(a+3)->GetBoundaryType()=='1') //24.02.2009. WW
					 elem->GetNode(a)->SetBoundaryType('1');

				// create a new tetrahedron
				Mesh_Group::CElem* new_elem = new Mesh_Group::CElem();
				new_elem->SetMark(true);
				new_elem->SetElementType(MshElemType::TETRAHEDRON);
				new_elem->SetPatchIndex(elem->GetPatchIndex());
				new_elem->SetBoundaryType('I');
				new_elem->SetNodesNumber(4);

				Math_Group::vec<Mesh_Group::CNode*> nodes(4);
				nodes[0] = mesh->nod_vector[elem->nodes_index[a]];
				nodes[1] = mesh->nod_vector[elem->nodes_index[b+3]];
				nodes[2] = mesh->nod_vector[elem->nodes_index[c+3]];
				nodes[3] = mesh->nod_vector[elem->nodes_index[c]];
				new_elem->SetNodes(nodes, true);

				new_elem->nodes_index.resize(4);
				for (size_t k=0; k<4; k++)
					new_elem->nodes_index[k] = elem->GetNode(k)->GetIndex();
				new_elems.push_back(new_elem);

				// change prism-element to 2nd tetrahedron
				elem->SetMark(true);
				elem->SetElementType(MshElemType::TETRAHEDRON);
				elem->SetNodesNumber(4);

				nodes[0] = mesh->nod_vector[elem->nodes_index[b]];
				nodes[1] = mesh->nod_vector[elem->nodes_index[a]];
				nodes[2] = mesh->nod_vector[elem->nodes_index[c]];
				nodes[3] = mesh->nod_vector[elem->nodes_index[b+3]];
				elem->SetNodes(nodes, true);

				elem->nodes_index.resize(4);
				for (size_t k=0; k<4; k++)
					elem->nodes_index[k] = elem->GetNode(k)->GetIndex();
				break;
			}
			case 2: // two nodes of the prism are marked false, i.e. create a tetrahedron element from the remaining 4 prism nodes
			{
				size_t a = false_node_idx[0];
				size_t b = (a+2)%3;
				size_t c = (a+1)%3;
				if (false_node_idx[1]==b)		a = c;
				else if(false_node_idx[1]==c)	a = b;

				elem->SetMark(true);
				elem->SetElementType(MshElemType::TETRAHEDRON);
				elem->SetNodesNumber(4);

				Math_Group::vec<Mesh_Group::CNode*> nodes(4);
				nodes[0] = mesh->nod_vector[elem->nodes_index[b]];
				nodes[1] = mesh->nod_vector[elem->nodes_index[a]];
				nodes[2] = mesh->nod_vector[elem->nodes_index[c]];
				nodes[3] = mesh->nod_vector[elem->nodes_index[a+3]];
				elem->SetNodes(nodes, true);

				elem->nodes_index.resize(4);
				for (size_t k=0; k<4; k++)
					elem->nodes_index[k] = elem->GetNode(k)->GetIndex();
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
	for(size_t i=0; i<new_elems.size(); i++)
		mesh->ele_vector.push_back(new_elems[i]);


	// correct indeces of elements and delete false elements
	std::vector<Mesh_Group::CElem*>::iterator beg_e = mesh->ele_vector.begin( ), last_e;
	long counter = 0;
	while ( beg_e != mesh->ele_vector.end() )
	{
		last_e = beg_e++;
		Mesh_Group::CElem* elem = *last_e;
		if (elem->GetMark())
		{
			elem->SetIndex(counter);
			counter++;
			for (int j=0; j<elem->GetVertexNumber(); j++)
			{
				if (!elem->GetNode(j)->GetMark())
				{
					Mesh_Group::CNode* node = elem->GetNode(j);
					node = mesh->nod_vector[elem->GetNode(j)->connected_nodes[0]];
				}
			}
		}
		else
		{
			delete elem;
			beg_e = mesh->ele_vector.erase(last_e);
		}
	}

	// correct indeces of nodes and delete false nodes
	counter = 0;
	std::vector<Mesh_Group::CNode*>::iterator beg = mesh->nod_vector.begin( ), last;
	while ( beg != mesh->nod_vector.end( ) )
	{
		last = beg++;
		Mesh_Group::CNode* node = *last;
		if (node->GetMark())
		{
			node->SetIndex(counter);
			node->connected_elements.clear();
			node->connected_nodes.clear();
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
	for(size_t i=0; i<nElems; i++)
	{
		for(int j=0; j< mesh->ele_vector[i]->GetVertexNumber(); j++)
			mesh->ele_vector[i]->nodes_index[j] = mesh->ele_vector[i]->GetNode(j)->GetIndex();
	}

	// delete elements if two of its nodes are identical (can this actually happen?????)
	beg_e = mesh->ele_vector.begin();
	counter = 0;
	bool flatf = false;
	while ( beg_e != mesh->ele_vector.end( ) )
	{
		last_e = beg_e++;
		Mesh_Group::CElem* elem = *last_e;

		//10.02.2009. WW !!!!!!!!!!!!!!!!!!!!!!
		for (int j=0; j<elem->GetVertexNumber(); j++)
		{
			flatf = false;
			for (int k=j; k<elem->GetVertexNumber(); k++)
			{
				if (elem->GetNodeIndex(j)==elem->GetNodeIndex(k))
				{
					flatf = true;
					break;
				}
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
		   for (int j=0; j<elem->GetVertexNumber(); j++)
		   {
			  if (!elem->GetNode(j)->GetMark())
			  {
				  Mesh_Group::CNode* node = elem->GetNode(j);
				  node = mesh->nod_vector[elem->GetNode(j)->connected_nodes[0]];
			  }
		   }
		}
	}

	mesh->ConnectedElements2Node();

	// delete nodes that are not connected to any element (can this happen???)
	counter = 0;
	beg = mesh->nod_vector.begin( );
	while ( beg != mesh->nod_vector.end( ) )
	{
		last = beg++;
		Mesh_Group::CNode* node = *last;
		if ( node->connected_elements.empty() )
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
	for (size_t i=0; i<nElems; i++)
	{
		for (int j=0; j<mesh->ele_vector[i]->GetVertexNumber(); j++)
			mesh->ele_vector[i]->nodes_index[j] = mesh->ele_vector[i]->GetNode(j)->GetIndex();
	}

	mesh->ConstructGrid();
}
