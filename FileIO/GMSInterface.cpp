/**
 * \file GMSInterface.cpp
 * 08/06/2010 KR Initial implementation
 *
 */

#include "GMSInterface.h"
#include "Station.h"
#include "Mesh.h"
#include "Node.h"
#include "Elements/Tet.h"
#include "Elements/Pyramid.h"
#include "Elements/Prism.h"
#include "MSHEnums.h"
#include "StringTools.h"
#include <fstream>

int GMSInterface::readBoreholesFromGMS(std::vector<GeoLib::Point*>* boreholes,
                                       const std::string &filename)
{
	double depth(-9999.0);
	std::string line(""), cName(""), sName("");
	std::list<std::string>::const_iterator it;
	GeoLib::Point* pnt = new GeoLib::Point();
	GeoLib::StationBorehole* newBorehole = NULL;
	std::ifstream in( filename.c_str() );

	if (!in.is_open())
	{
		std::cout << "GMSInterface::readBoreholeFromGMS() - Could not open file...\n";
		return 0;
	}

	/* skipping first line because it contains field names */
	getline(in, line);

	/* read all stations */
	while ( getline(in, line) )
	{
		std::list<std::string> fields = splitString(line, '\t');

		if (fields.size() >= 5)
		{
			if (fields.begin()->compare(cName) == 0) // add new layer
			{
				it = fields.begin();
				(*pnt)[0] = strtod((++it)->c_str(), 0);
				(*pnt)[1] = strtod((++it)->c_str(), 0);
				(*pnt)[2] = strtod((++it)->c_str(), 0);

				// check if current layer has a thickness of 0.0.
				// if so skip it since it will mess with the vtk-visualisation later on!
				if ((*pnt)[2] != depth)
				{
					newBorehole->addSoilLayer((*pnt)[0], (*pnt)[1], (*pnt)[2], sName);
					sName = (*(++it));
					depth = (*pnt)[2];
				}
				else
					std::cout << "Warning: Skipped layer \"" << sName << "\" in borehole \"" 
					          << cName << "\" because of thickness 0.0." << std::endl;
			}
			else // add new borehole
			{
				if (newBorehole != NULL)
				{
					newBorehole->setDepth((*newBorehole)[2] - depth);
					boreholes->push_back(newBorehole);
				}
				cName = *fields.begin();
				it = fields.begin();
				(*pnt)[0] = strtod((++it)->c_str(), 0);
				(*pnt)[1] = strtod((++it)->c_str(), 0);
				(*pnt)[2] = strtod((++it)->c_str(), 0);
				sName = (*(++it));
				newBorehole = GeoLib::StationBorehole::createStation(cName, (*pnt)[0], (*pnt)[1], (*pnt)[2], 0);
				depth = (*pnt)[2];
			}
		}
		else
			std::cout <<
			"GMSInterface::readBoreholeFromGMS() - Error reading format..." <<
			std::endl;
	}
	// write the last borehole from the file
	if (newBorehole != NULL)
	{
		newBorehole->setDepth((*newBorehole)[2] - depth);
		boreholes->push_back(newBorehole);
	}
	in.close();

	if (boreholes->empty())
		return 0;
	return 1;
}

/*
   // all boreholes to GMS which each borehole in a single file
   void StationIO::writeBoreholesToGMS(const std::vector<GeoLib::Point*> *stations)
   {
    //std::vector<std::string> soilID(1);
    std::vector<std::string> soilID = readSoilIDfromFile("d:/BodeTimeline.txt");
    for (size_t i=0; i<stations->size(); i++)
        StationIO::writeBoreholeToGMS(static_cast<GeoLib::StationBorehole*>((*stations)[i]), std::string("Borehole-" + static_cast<GeoLib::StationBorehole*>((*stations)[i])->getName() + ".txt"), soilID);
    StationIO::writeSoilIDTable(soilID, "SoilIDReference.txt");
   }
 */
void GMSInterface::writeBoreholesToGMS(const std::vector<GeoLib::Point*>* stations,
                                       const std::string &filename)
{
	std::ofstream out( filename.c_str(), std::ios::out );
	size_t idx = 0;
	std::vector<std::string> soilID = readSoilIDfromFile("d:/BodeTimeline.txt");

	// write header
	out << "name" << "\t" << std::fixed << "X" << "\t" << "Y"  << "\t" << "Z" <<  "\t" <<
	"soilID" << std::endl;

	for (size_t j = 0; j < stations->size(); j++)
	{
		GeoLib::StationBorehole* station =
		        static_cast<GeoLib::StationBorehole*>((*stations)[j]);
		std::vector<GeoLib::Point*> profile = station->getProfile();
		std::vector<std::string> soilNames  = station->getSoilNames();

		size_t nLayers = profile.size();
		for (size_t i = 1; i < nLayers; i++)
		{
			if ( (i > 1) && (soilNames[i].compare(soilNames[i - 1]) == 0) )
				continue;
			idx = getSoilID(soilID, soilNames[i]);

			out << station->getName() << "\t" << std::fixed <<
			(*(profile[i - 1]))[0] << "\t"
			    << (*(profile[i - 1]))[1]  << "\t" << (*(profile[i - 1]))[2] <<  "\t"
			    << idx << std::endl;
		}
		out << station->getName() << "\t" << std::fixed <<
		(*(profile[nLayers - 1]))[0] << "\t"
		    << (*(profile[nLayers -
		              1]))[1]  << "\t" << (*(profile[nLayers - 1]))[2] <<  "\t"
		    << idx << std::endl; // this line marks the end of the borehole
	}

	out.close();
	GMSInterface::writeSoilIDTable(soilID, "d:/SoilIDReference.txt");
}

int GMSInterface::writeBoreholeToGMS(const GeoLib::StationBorehole* station,
                                     const std::string &filename,
                                     std::vector<std::string> &soilID)
{
	std::ofstream out( filename.c_str(), std::ios::out );
	size_t idx = 0;

	// write header
	out << "name" << "\t" << std::fixed << "X" << "\t" << "Y"  << "\t" << "Z" <<  "\t" <<
	"soilID" << std::endl;

	std::vector<GeoLib::Point*> profile = station->getProfile();
	std::vector<std::string> soilNames  = station->getSoilNames();

	// write table
	size_t nLayers = profile.size();
	for (size_t i = 1; i < nLayers; i++)
	{
		if ( (i > 1) && (soilNames[i].compare(soilNames[i - 1]) == 0) )
			continue;
		idx = getSoilID(soilID, soilNames[i]);

		out << station->getName() << "\t" << std::fixed << (*(profile[i - 1]))[0] << "\t"
		    << (*(profile[i - 1]))[1]  << "\t" << (*(profile[i - 1]))[2] <<  "\t"
		    << idx << std::endl;
	}
	out << station->getName() << "\t" << std::fixed << (*(profile[nLayers - 1]))[0] << "\t"
	    << (*(profile[nLayers - 1]))[1]  << "\t" << (*(profile[nLayers - 1]))[2] <<  "\t"
	    << idx << std::endl;        // this line marks the end of the borehole
	out.close();

	return 1;
}

size_t GMSInterface::getSoilID(std::vector<std::string> &soilID, std::string &soilName)
{
	for (size_t j = 0; j < soilID.size(); j++)
		if (soilID[j].compare(soilName) == 0)
			return j;
	soilID.push_back(soilName);
	return soilID.size() - 1;
}

int GMSInterface::writeSoilIDTable(const std::vector<std::string> &soilID,
                                   const std::string &filename)
{
	std::ofstream out( filename.c_str(), std::ios::out );

	// write header
	out << "ID" << "\t" << std::fixed << "Soil name" << std::endl;

	// write table
	size_t nIDs = soilID.size();
	for (size_t i = 0; i < nIDs; i++)
		out << i << "\t" << std::fixed << soilID[i] << "\t" << std::endl;
	out.close();

	return 1;
}

std::vector<std::string> GMSInterface::readSoilIDfromFile(const std::string &filename)
{
	std::vector<std::string> soilID;
	std::string line;

	std::ifstream in( filename.c_str() );

	if (in.is_open())
		while ( getline(in, line) )
		{
			trim(line);
			soilID.push_back(line);
		}
	in.close();

	return soilID;
}

MeshLib::Mesh* GMSInterface::readGMS3DMMesh(std::string filename)
{
	std::string line("");

	std::ifstream in(filename.c_str());
	if (!in.is_open())
	{
		std::cout << "GMSInterface::readGMS3DMMesh() - Could not open file..." << std::endl;
		return NULL;
	}

	// Read data from file
	getline(in, line); // "MESH3D"
	if (line.compare("MESH3D") != 0)
	{
		std::cout << "GMSInterface::readGMS3DMMesh() - Could not read expected file header..." << std::endl;
		return NULL;
	}

	std::cout << "Read GMS-3DM mesh ... ";
	std::vector<MeshLib::Node*> nodes;
	std::vector<MeshLib::Element*> elements;
	std::map<unsigned, unsigned> id_map;

	// elements are listed before nodes in 3dm-format, therefore
	// traverse file twice and read first nodes and then elements
	std::string dummy;
	unsigned id(0), count(0);
	double x[3];
	// read nodes
	while ( getline(in, line) )
	{
		if (line[0] == 'N') // "ND" for Node
		{
			std::stringstream str(line);
			str >> dummy >> id >> x[0] >> x[1] >> x[2];
			MeshLib::Node* node = new MeshLib::Node(x, id);
			id_map.insert(std::pair<unsigned, unsigned>(id,count++));
			nodes.push_back(node);
		}
	}
	in.close();

	// NOTE: Element types E8H (Hex), E4Q (Quad), E3T (Tri) are not implemented yet
	// read elements
	in.open(filename.c_str());
	getline(in, line); // "MESH3D"
	unsigned node_idx[6], mat_id;
	while ( getline(in, line) )
	{
		MeshLib::Element* elem (NULL);
		std::string element_id(line.substr(0,3));
		std::stringstream str(line);

		if (element_id.compare("E6W") == 0)	// Prism
		{
			str >> dummy >> id >> node_idx[0] >> node_idx[1] >> node_idx[2] >> node_idx[3] 
			    >> node_idx[4] >> node_idx[5] >> mat_id;
			elem = new MeshLib::Prism(nodes[id_map.find(node_idx[0])->second], nodes[id_map.find(node_idx[1])->second], 
									  nodes[id_map.find(node_idx[2])->second], nodes[id_map.find(node_idx[3])->second], 
									  nodes[id_map.find(node_idx[4])->second], nodes[id_map.find(node_idx[5])->second], mat_id);
			elements.push_back(elem);
		}
		else if (element_id.compare("E4T") == 0) // Tet
		{
			str >> dummy >> id >> node_idx[0] >> node_idx[1] >> node_idx[2] >> node_idx[3] >> mat_id;
			elem = new MeshLib::Tet(nodes[id_map.find(node_idx[0])->second], nodes[id_map.find(node_idx[1])->second], 
				                    nodes[id_map.find(node_idx[2])->second], nodes[id_map.find(node_idx[3])->second], mat_id);
			elements.push_back(elem);
		}
		else if ((element_id.compare("E4P") == 0) || (element_id.compare("E5P") == 0)) // Pyramid (both do exist for some reason)
		{
			str >> dummy >> id >> node_idx[0] >> node_idx[1] >> node_idx[2] >> node_idx[3] >> node_idx[4] >> mat_id;
			elem = new MeshLib::Pyramid(nodes[id_map.find(node_idx[0])->second], nodes[id_map.find(node_idx[1])->second], 
				                        nodes[id_map.find(node_idx[2])->second], nodes[id_map.find(node_idx[3])->second], 
										nodes[id_map.find(node_idx[4])->second], mat_id);
			elements.push_back(elem);
		}
		else if (element_id.compare("ND ") == 0) // Node
		{
			continue;	// skip because nodes have already been read
		}
		else //default
		{
			std::cout << std::endl << "GMSInterface::readGMS3DMMesh() - Element type \"" << element_id << "\"not recognised ..." << std::endl;
			return NULL;
		}
	}

	in.close();
	std::cout << "finished" << std::endl;

	std::string mesh_name (BaseLib::getFileNameFromPath(filename));
	return new MeshLib::Mesh(mesh_name, nodes, elements);
}
