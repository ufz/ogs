/**
 * \file removeMeshNodes.cpp
 * 2012/03/07 KR Initial implementation
 */

#include <QApplication>
#include "logog/include/logog.hpp"
#include "LogogSimpleFormatter.h"
#include "readMeshFromFile.h"
#include "Legacy/MeshIO.h"
#include "AABB.h"
#include "Mesh.h"
#include "Node.h"
#include "MeshEditing/removeMeshNodes.h"
#include "MathTools.h"

int find_closest_point(MeshLib::Node const*const point, std::vector<MeshLib::Node*> const& nodes, double const& max_dist)
{
	const std::size_t nNodes (nodes.size());
	double sqr_shortest_dist (max_dist*2);
	int idx = (sqr_shortest_dist<max_dist) ? 0 : -1;
	const MeshLib::Node p (*point);

	for (unsigned i=0; i<nNodes; i++) 
	{
		double sqr_dist ((p[0]-(*nodes[i])[0])*(p[0]-(*nodes[i])[0]));
		if (sqr_dist < max_dist)
		{
			sqr_dist += ((p[1]-(*nodes[i])[1])*(p[1]-(*nodes[i])[1]));
			if (sqr_dist < max_dist && sqr_dist < sqr_shortest_dist)
			{
				sqr_shortest_dist = sqr_dist;
				idx = i;
			}
		}
	}

	return idx;
}

bool containsPoint(MeshLib::Node const& pnt, MeshLib::Node const& min, MeshLib::Node const& max) 
{
	if (pnt[0] < min[0] || max[0] < pnt[0]) return false;
	if (pnt[1] < min[1] || max[1] < pnt[1]) return false;
	return true;
}

int main (int argc, char* argv[])
{
	QApplication app(argc, argv, false);
	LOGOG_INITIALIZE();
	logog::Cout* logogCout = new logog::Cout;
	BaseLib::LogogSimpleFormatter* formatter = new BaseLib::LogogSimpleFormatter;
	logogCout->SetFormatter(*formatter);

	std::vector<std::string> keywords;
	keywords.push_back("-ALL");
	keywords.push_back("-MESH");
	keywords.push_back("-LOWPASS");

	if (argc < 3)
	{
		std::cout << "Moves mesh nodes and connected elements either by a given value or based on a list." << std::endl;
		std::cout << std::endl;
		std::cout << "Usage: " << argv[0] << " <msh-file.msh> <keyword> [<value1>] [<value2>]" << std::endl;
		std::cout << "Available keywords:" << std::endl;
		//for (size_t i=0; i<keywords.size(); i++)
		std::cout << "\t" << "-ALL <value1> <value2> : changes the elevation of all mesh nodes by <value2> in direction <value1> [x,y,z]." << std::endl;
		std::cout << "\t" << "-MESH <value1> <value2> : changes the elevation of mesh nodes based on a second mesh <value1> with a search range of <value2>." << std::endl;
		std::cout << "\t" << "-LOWPASS : applies a lowpass filter over node elevation using directly connected nodes." << std::endl;
		return -1;
	}
	
	const std::string msh_name(argv[1]);
	const std::string current_key(argv[2]);
	//const std::string msh_name("D:\\rappbode-2013-03-03--30m_lowpass_new_new.msh");
	//const std::string current_key("-MESH");

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
	//std::vector<size_t> del_nodes;

	// Start keyword-specific selection of nodes

	// moves the elevation of all nodes by value
	if (current_key.compare("-ALL")==0)
	{
		if (argc < 5)
		{
			std::cout << "Missing parameter..." << std::endl;
			return -1;
		}
		const std::string dir(argv[3]);
		unsigned idx = (dir.compare("x") == 0) ? 0 : (dir.compare("y") == 0) ? 1 : 2;
		const double value(strtod(argv[4],0));
		std::cout << "Moving all mesh nodes by " << value << " in direction " << idx << " (" << dir << ")..." << std::endl;
		//double value(-10);
		const size_t nNodes(mesh->getNNodes());
		std::vector<MeshLib::Node*> nodes (mesh->getNodes());
		for (size_t i=0; i<nNodes; i++)
		{
			(*nodes[i])[idx] += value;
		}
	}

	// maps the elevation of mesh nodes according to a ground truth mesh whenever nodes exist within max_dist
	if (current_key.compare("-MESH")==0)
	{
		if (argc < 5)
		{
			std::cout << "Missing parameter..." << std::endl;
			return -1;
		}
		const std::string value (argv[3]);
		double max_dist(pow(strtod(argv[4],0), 2));
		//const std::string value("D:\\Rappbodevorsperre_elevation440m.msh");
		//double max_dist (25.0);	// squared maximum distance at which reference points are used
		double offset (0.0); // additional offset for elevation (should be 0)
		MeshLib::Mesh* ground_truth (FileIO::readMeshFromFile(value));
		const std::vector<MeshLib::Node*> ground_truth_nodes (ground_truth->getNodes());
		GeoLib::AABB<MeshLib::Node> bounding_box(ground_truth_nodes.begin(), ground_truth_nodes.end());
		const MeshLib::Node min (bounding_box.getMinPoint());
		const MeshLib::Node max (bounding_box.getMaxPoint());
		
		const size_t nNodes(mesh->getNNodes());
		std::vector<MeshLib::Node*> nodes (mesh->getNodes());

		for (size_t i=0; i<nNodes; i++)
		{
			bool is_inside (containsPoint(*nodes[i], min, max));
			if (is_inside)
			{
				int idx = find_closest_point(nodes[i], ground_truth_nodes, max_dist);
				if (idx>=0)
					(*nodes[i])[2] = (*(ground_truth_nodes[idx]))[2]-offset;
			}
		}
	}

	// a simple lowpass filter for the elevation of mesh nodes using the elevation of each node 
	// weighted by 2 and the elevation of each connected node weighted by 1
	if (current_key.compare("-LOWPASS")==0)
	{
		const size_t nNodes(mesh->getNNodes());
		std::vector<MeshLib::Node*> nodes (mesh->getNodes());

		std::vector<double> elevation(nNodes);
		for (size_t i=0; i<nNodes; i++)
			elevation[i] = (*nodes[i])[2];

		for (size_t i=0; i<nNodes; i++)
		{
			const std::vector<MeshLib::Node*> conn_nodes (nodes[i]->getConnectedNodes());
			const unsigned nConnNodes (conn_nodes.size());
			elevation[i] = (2*(*nodes[i])[2]);
			for (size_t j=0; j<nConnNodes; ++j)
				elevation[i] += (*conn_nodes[j])[2];
			elevation[i] /= (nConnNodes+2);
		}

		for (size_t i=0; i<nNodes; i++)
			(*nodes[i])[2] = elevation[i];
	}
	/**** add other keywords here ****/

	FileIO::Legacy::MeshIO meshIO;
	meshIO.setMesh(mesh);
	meshIO.setPrecision(9);
	meshIO.writeToFile(msh_name.substr(0, msh_name.length()-4) + "_new.msh");
	delete mesh;
	return 1;

}



