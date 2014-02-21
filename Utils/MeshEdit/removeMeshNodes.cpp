/**
 * \file removeMeshNodes.cpp
 * 2012/03/07 KR Initial implementation
 */

#include <QApplication>
#include "readMeshFromFile.h"
#include "Legacy/MeshIO.h"
#include "Mesh.h"
#include "Node.h"
#include "MeshEditing/removeMeshNodes.h"

int main (int argc, char* argv[])
{
	QApplication app(argc, argv, false);

	std::vector<std::string> keywords;
	keywords.push_back("-Z_ABOVE");
	keywords.push_back("-Z_BELOW");

	if (argc != 4)
	{
		std::cout << "Removes mesh nodes and connected elements based on the given criterium." << std::endl;
		std::cout << std::endl;
		std::cout << "Usage: " << argv[0] << " <msh-file.msh> <keyword> <value>" << std::endl;
		std::cout << "Available keywords:" << std::endl;
		for (size_t i=0; i<keywords.size(); i++)
			std::cout << "\t" << keywords[i] << std::endl;
		return -1;
	}

	const std::string msh_name(argv[1]);
	const std::string current_key(argv[2]);
	double value(strtod(argv[3],0));

	if (msh_name.substr(msh_name.length()-4, 4).compare(".msh") != 0)
	{
		std::cout << "Error: Parameter 1 should be a msh-file" << std::endl;
		std::cout << "Usage: " << argv[0] << " <msh-file.gml> <keyword> <value>" << std::endl;
		return -1;
	}

	bool is_keyword(false);
	for (size_t i=0; i<keywords.size(); i++)
		if (current_key.compare(keywords[i])==0)
		{
			is_keyword = true;
			break;
		}

	if (!is_keyword)
	{
		std::cout << "Keyword not recognised. Available keywords:" << std::endl;
		for (size_t i=0; i<keywords.size(); i++)
			std::cout << keywords[i] << std::endl;
		return -1;
	}

	MeshLib::Mesh* mesh (FileIO::readMeshFromFile(msh_name));
	std::vector<size_t> del_nodes;

	// Start keyword-specific selection of nodes
	if ((current_key.compare("-Z_BELOW")==0) || (current_key.compare("-Z_ABOVE")==0))
	{
		int factor(1);
		if (current_key.compare("-Z_ABOVE")==0)
		{
			factor=-1;
			value*=-1;
		}

		const size_t nNodes(mesh->getNNodes());
		for (size_t i=0; i<nNodes; i++)
		{
			const MeshLib::Node* node = mesh->getNode(i);
			const double* coords(node->getCoords());
			if ((factor*coords[2]) < value)
				del_nodes.push_back(i);
		}
	}

	/**** add other keywords here ****/

	// remove nodes and write new file
	MeshLib::removeMeshNodes(*mesh, del_nodes);
	
	FileIO::Legacy::MeshIO meshIO;
	meshIO.setMesh(mesh);
	meshIO.setPrecision(9);
	meshIO.writeToFile(msh_name.substr(0, msh_name.length()-4) + "_new.msh");
	delete mesh;
	return 1;
}



