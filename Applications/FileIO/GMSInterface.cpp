/**
 * \file
 * \author Karsten Rink
 * \date   2010-06-08
 * \brief  Implementation of the GMSInterface class.
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * @file GMSInterface.cpp
 * @date 2010-06-08
 * @author Karsten Rink
 */

#include "GMSInterface.h"

#include <fstream>
#include <logog/include/logog.hpp>

#include "BaseLib/FileTools.h"
#include "BaseLib/StringTools.h"

#include "GeoLib/StationBorehole.h"

#include "MeshLib/Elements/Prism.h"
#include "MeshLib/Elements/Pyramid.h"
#include "MeshLib/Elements/Tet.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshEnums.h"
#include "MeshLib/Node.h"

namespace FileIO
{

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
        ERR("GMSInterface::readBoreholeFromGMS(): Could not open file %s.", filename.c_str());
        return 0;
    }

    /* skipping first line because it contains field names */
    getline(in, line);

    /* read all stations */
    while ( getline(in, line) )
    {
        std::list<std::string> fields = BaseLib::splitString(line, '\t');

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
                    newBorehole->addSoilLayer((*pnt)[0],
                                              (*pnt)[1],
                                              (*pnt)[2],
                                              sName);
                    sName = (*(++it));
                    depth = (*pnt)[2];
                }
                else
                    WARN("GMSInterface::readBoreholeFromGMS(): Skipped layer \"%s\" in borehole \"%s\" because of thickness 0.0.",
                         sName.c_str(), cName.c_str());
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
                newBorehole =
                        GeoLib::StationBorehole::createStation(cName, (*pnt)[0],
                                                               (*pnt)[1],
                                                               (*pnt)[2], 0);
                depth = (*pnt)[2];
            }
        }
        else
            ERR("GMSInterface::readBoreholeFromGMS(): Error reading format.");
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
    for (std::size_t i=0; i<stations->size(); i++)
        StationIO::writeBoreholeToGMS(static_cast<GeoLib::StationBorehole*>((*stations)[i]), std::string("Borehole-" + static_cast<GeoLib::StationBorehole*>((*stations)[i])->getName() + ".txt"), soilID);
    StationIO::writeSoilIDTable(soilID, "SoilIDReference.txt");
   }
 */
void GMSInterface::writeBoreholesToGMS(const std::vector<GeoLib::Point*>* stations,
                                       const std::string &filename)
{
    std::ofstream out( filename.c_str(), std::ios::out );
    //std::vector<std::string> soilID = readSoilIDfromFile("d:/BodeTimeline.txt");

    // write header
    out << "name" << "\t" << std::fixed << "X" << "\t" << "Y"  << "\t" << "Z" <<  "\t" <<
    "soilID" << "\n";

    for (auto station_as_point : *stations)
    {
        GeoLib::StationBorehole* station =
                static_cast<GeoLib::StationBorehole*>(station_as_point);
        std::vector<GeoLib::Point*> profile = station->getProfile();
        std::vector<std::string> soilNames  = station->getSoilNames();
        //std::size_t idx = 0;
        std::string current_soil_name("");

        std::size_t nLayers = profile.size();
        for (std::size_t i = 1; i < nLayers; i++)
        {
            if ( (i > 1) && (soilNames[i].compare(soilNames[i - 1]) == 0) )
                continue;
            //idx = getSoilID(soilID, soilNames[i]);
            current_soil_name = soilNames[i];

            out << station->getName() << "\t" << std::fixed
                << (*(profile[i - 1]))[0] << "\t"
                << (*(profile[i - 1]))[1]  << "\t"
                << (*(profile[i - 1]))[2] <<  "\t"
                << current_soil_name/*idx*/ << "\n";
        }
        out << station->getName() << "\t" << std::fixed <<
        (*(profile[nLayers - 1]))[0] << "\t"
            << (*(profile[nLayers - 1]))[1]  << "\t"
            << (*(profile[nLayers - 1]))[2] <<  "\t"
            << current_soil_name << "\n"; // this line marks the end of the borehole
    }

    out.close();
    //GMSInterface::writeSoilIDTable(soilID, "d:/SoilIDReference.txt");
}

std::size_t GMSInterface::getSoilID(std::vector<std::string> &soilID, std::string &soilName)
{
    for (std::size_t j = 0; j < soilID.size(); j++)
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
    out << "ID" << "\t" << std::fixed << "Soil name" << "\n";

    // write table
    std::size_t nIDs = soilID.size();
    for (std::size_t i = 0; i < nIDs; i++)
        out << i << "\t" << std::fixed << soilID[i] << "\t" << "\n";
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
            BaseLib::trim(line);
            soilID.push_back(line);
        }
    in.close();

    return soilID;
}

MeshLib::Mesh* GMSInterface::readGMS3DMMesh(const std::string &filename)
{
    std::string line("");

    std::ifstream in(filename.c_str());
    if (!in.is_open())
    {
        ERR("GMSInterface::readGMS3DMMesh(): Could not open file %s.", filename.c_str());
        return NULL;
    }

    // Read data from file
    getline(in, line); // "MESH3D"
    if (line.compare("MESH3D") != 0)
    {
        ERR("GMSInterface::readGMS3DMMesh(): Could not read expected file header.");
        return NULL;
    }

    INFO("GMSInterface::readGMS3DMMesh(): Read GMS-3DM mesh.");
    std::vector<MeshLib::Node*> nodes;
    std::vector<MeshLib::Element*> elements;
    std::vector<int> mat_ids;
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
    unsigned node_idx[6];
    int mat_id(0);
    while ( getline(in, line) )
    {
        std::string element_id(line.substr(0,3));
        std::stringstream str(line);

        if (element_id.compare("E6W") == 0) // Prism
        {
            str >> dummy >> id >> node_idx[0] >> node_idx[1] >> node_idx[2] >>
            node_idx[3]
            >> node_idx[4] >> node_idx[5] >> mat_id;
            MeshLib::Node** prism_nodes = new MeshLib::Node*[6];
            for (unsigned k(0); k<6; k++) {
                prism_nodes[k] = nodes[id_map.find(node_idx[k])->second];
            }
            elements.push_back(new MeshLib::Prism(prism_nodes));
            mat_ids.push_back(mat_id);
        }
        else if (element_id.compare("E4T") == 0) // Tet
        {
            str >> dummy >> id >> node_idx[0] >> node_idx[1] >> node_idx[2] >>
            node_idx[3] >> mat_id;
            MeshLib::Node** tet_nodes = new MeshLib::Node*[4];
            for (unsigned k(0); k<4; k++) {
                tet_nodes[k] = nodes[id_map.find(node_idx[k])->second];
            }
            elements.push_back(new MeshLib::Tet(tet_nodes));
            mat_ids.push_back(mat_id);
        }
        else if ((element_id.compare("E4P") == 0) || (element_id.compare("E5P") == 0)) // Pyramid (both do exist for some reason)
        {
            str >> dummy >> id >> node_idx[0] >> node_idx[1] >> node_idx[2] >>
            node_idx[3] >> node_idx[4] >> mat_id;
            MeshLib::Node** pyramid_nodes = new MeshLib::Node*[5];
            for (unsigned k(0); k<5; k++) {
                pyramid_nodes[k] = nodes[id_map.find(node_idx[k])->second];
            }
            elements.push_back(new MeshLib::Pyramid(pyramid_nodes));
            mat_ids.push_back(mat_id);
        }
        else if (element_id.compare("ND ") == 0) // Node

            continue; // skip because nodes have already been read
        else //default
        {
            WARN("GMSInterface::readGMS3DMMesh() - Element type \"%s\" not recognised.",
                 element_id.c_str());
            return NULL;
        }
    }

    in.close();
    INFO("GMSInterface::readGMS3DMMesh(): finished.");

    const std::string mesh_name = BaseLib::extractBaseNameWithoutExtension(filename);
    MeshLib::Properties properties;
    if (mat_ids.size() == elements.size())
    {
        auto* const opt_pv = properties.createNewPropertyVector<int>(
            "MaterialIDs", MeshLib::MeshItemType::Cell);
        if (!opt_pv)
        {
            ERR("Could not create PropertyVector for material ids.");
            for (auto element : elements)
                delete element;
            for (auto node : nodes)
                delete node;
            return nullptr;
        }
        opt_pv->reserve(mat_ids.size());
        std::copy(mat_ids.cbegin(), mat_ids.cend(),
                  std::back_inserter(*opt_pv));
    }
    else
        ERR ("Ignoring Material IDs information (does not match number of elements).");
    return new MeshLib::Mesh(mesh_name, nodes, elements, properties);
}

}
