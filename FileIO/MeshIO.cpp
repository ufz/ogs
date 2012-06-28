/**
 * \file MeshIO.cpp
 *
 * Created on 2012-05-08 by Karsten Rink
 */

#include "GEOObjects.h"
#include "MeshIO.h"
#include "Node.h"
#include "Elements/Edge.h"
#include "Elements/Tri.h"
#include "Elements/Quad.h"
#include "Elements/Tet.h"
#include "Elements/Hex.h"
#include "Elements/Pyramid.h"
#include "Elements/Prism.h"

#include <iomanip>
#include <sstream>

namespace FileIO
{

MeshIO::MeshIO()
: _mesh(NULL)
{
}

MeshLib::Mesh* MeshIO::loadMeshFromFile(const std::string& file_name)
{
	std::cout << "Read mesh ... " << std::endl;

	std::ifstream in (file_name.c_str(),std::ios::in);
	if (!in.is_open())
	{
		std::cout << "CFEMesh::FEMRead() - Could not open file...\n";
		return NULL;
	}

	std::string line_string ("");
	getline(in, line_string);

	std::vector<MeshLib::Node*> nodes;
	std::vector<MeshLib::Element*> elements;

	if(line_string.find("#FEM_MSH") != std::string::npos) // OGS mesh file
	{
		double edge_length[2] = { std::numeric_limits<double>::max(), std::numeric_limits<double>::min() };
		while (!in.eof())
		{
			getline(in, line_string);

			// check keywords
			if (line_string.find("#STOP") != std::string::npos)
				break;
			else if (line_string.find("$NODES") != std::string::npos)
			{
				double x, y, z, double_dummy;
				unsigned nNodes, idx;
				in >> nNodes >> std::ws;
				std::string s;
				std::ios::pos_type position = in.tellg();
				for (unsigned i = 0; i < nNodes; i++)
				{
					in >> idx >> x >> y >> z;
					MeshLib::Node* node(new MeshLib::Node(x, y, z, nodes.size()));
					nodes.push_back(node);
					position = in.tellg();
					in >> s;
					if (s.find("$AREA") != std::string::npos)
						in >> double_dummy;
					else
						in.seekg(position, std::ios::beg);
					in >> std::ws;
				}
			}
			else if (line_string.find("$ELEMENTS") != std::string::npos)
			{
				unsigned nElements;
				in >> nElements >> std::ws;
				for (unsigned i = 0; i < nElements; i++)
				{
					getline(in, line_string);

					size_t elem_idx (elements.size());
					elements.push_back(readElement(line_string, nodes));

					double elem_min_length, elem_max_length;
					elements[elem_idx]->computeSqrEdgeLengthRange(elem_min_length, elem_max_length);
					edge_length[0] = (elem_min_length<edge_length[0]) ? elem_min_length : edge_length[0];
					edge_length[1] = (elem_max_length>edge_length[1]) ? elem_max_length : edge_length[1];
				}
			}
		}

		MeshLib::Mesh* mesh (new MeshLib::Mesh(file_name, nodes, elements));
		mesh->setEdgeLengthRange(sqrt(edge_length[0]), sqrt(edge_length[1]));

		std::cout << "Nr. Nodes: " << nodes.size() << std::endl;
		std::cout << "Nr. Elements: " << elements.size() << std::endl;

		in.close();
		return mesh;
	}
	else
	{
		in.close();
		return NULL;
	}
}

MeshLib::Element* MeshIO::readElement(const std::string& line, const std::vector<MeshLib::Node*> &nodes)
{
	std::stringstream ss (line);
	std::string elem_type_str;
	unsigned index, patch_index;
	ss >> index >> patch_index >> elem_type_str;

	MshElemType::type elem_type (String2MshElemType(elem_type_str));
	unsigned* idx = new unsigned[8];

	MeshLib::Element* elem;

	switch(elem_type)
	{
	case MshElemType::LINE:
		for (int i = 0; i < 2; i++)
			ss >> idx[i];
		elem = new MeshLib::Edge(nodes[idx[0]], nodes[idx[1]], patch_index);
		break;
	case MshElemType::TRIANGLE:
		for (int i = 0; i < 3; i++)
			ss >> idx[i];
		elem = new MeshLib::Tri(nodes[idx[0]], nodes[idx[1]], nodes[idx[2]], patch_index);
		break;
	case MshElemType::QUAD:
		for (int i = 0; i < 4; i++)
			ss >> idx[i];
		elem = new MeshLib::Quad(nodes[idx[0]], nodes[idx[1]], nodes[idx[2]], nodes[idx[3]], patch_index);
		break;
	case MshElemType::TETRAHEDRON:
		for (int i = 0; i < 4; i++)
			ss >> idx[i];
		elem = new MeshLib::Tet(nodes[idx[0]], nodes[idx[1]], nodes[idx[2]], nodes[idx[3]], patch_index);
		break;
	case MshElemType::HEXAHEDRON:
		for (int i = 0; i < 8; i++)
			ss >> idx[i];
		elem = new MeshLib::Hex(nodes[idx[0]], nodes[idx[1]], nodes[idx[2]], nodes[idx[3]], nodes[idx[4]], nodes[idx[5]], nodes[idx[6]], nodes[idx[7]], patch_index);
		break;
	case MshElemType::PYRAMID:
		for (int i = 0; i < 5; i++)
			ss >> idx[i];
		elem = new MeshLib::Pyramid(nodes[idx[0]], nodes[idx[1]], nodes[idx[2]], nodes[idx[3]], nodes[idx[4]], patch_index);
		break;
	case MshElemType::PRISM:
		for (int i = 0; i < 6; i++)
			ss >> idx[i];
		elem = new MeshLib::Prism(nodes[idx[0]], nodes[idx[1]], nodes[idx[2]], nodes[idx[3]], nodes[idx[4]], nodes[idx[5]], patch_index);
		break;
	default:
		elem = NULL;
	}

	return elem;
}

int MeshIO::write(std::ostream &out)
{
	(void)out;
	/*
	if(!_mesh)
	{
		std::cout << "OGSMeshIO cannot write: no mesh set!" << std::endl;
		return 0;
	}

	setPrecision(9);

	out << "#FEM_MSH" << std::endl;

	out << "$PCS_TYPE" << std::endl << "  " << _mesh->pcs_name << std::endl;

	out << "$NODES" << std::endl << "  ";
	const size_t n_nodes(_mesh->GetNodesNumber(false));
	out << n_nodes << std::endl;
	for (size_t i(0); i < n_nodes; i++)
	{
		double const* const coords (_mesh->nod_vector[i]->getCoords());
		out << i << " " << coords[0] << " " << coords[1] << " " << coords[2] << std::endl;
	}

	out << "$ELEMENTS" << std::endl << "  ";
	writeElementsExceptLines(_mesh->ele_vector, out);

	out << " $LAYER" << std::endl;
	out << "  ";
	out << _mesh->_n_msh_layer << std::endl;
	out << "#STOP" << std::endl;
	*/
	return 1;
}

void MeshIO::setMesh(const MeshLib::Mesh* mesh)
{
	_mesh = mesh;
}

} // end namespace FileIO
