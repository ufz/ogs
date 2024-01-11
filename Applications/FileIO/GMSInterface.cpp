/**
 * \file
 * \author Karsten Rink
 * \date   2010-06-08
 * \brief  Implementation of the GMSInterface class.
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "GMSInterface.h"

#include <fstream>
#include <list>

#include "BaseLib/FileTools.h"
#include "BaseLib/Logging.h"
#include "BaseLib/StringTools.h"
#include "GeoLib/StationBorehole.h"
#include "MeshLib/Elements/Prism.h"
#include "MeshLib/Elements/Pyramid.h"
#include "MeshLib/Elements/Tet.h"
#include "MeshLib/Elements/Tri.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshEnums.h"
#include "MeshLib/Node.h"

namespace
{
template <typename It>
std::array<double, 3> parsePointCoordinates(It& it)
{
    return {std::strtod((++it)->c_str(), nullptr),
            std::strtod((++it)->c_str(), nullptr),
            std::strtod((++it)->c_str(), nullptr)};
}
}  // namespace

namespace FileIO
{
int GMSInterface::readBoreholesFromGMS(std::vector<GeoLib::Point*>& boreholes,
                                       const std::string& filename)
{
    std::ifstream in(filename.c_str());
    if (!in.is_open())
    {
        ERR("GMSInterface::readBoreholeFromGMS(): Could not open file {:s}.",
            filename);
        return 0;
    }

    double depth(-9999.0);
    std::string line;
    std::string cName;
    std::string sName;
    std::list<std::string>::const_iterator it;
    GeoLib::StationBorehole* newBorehole = nullptr;

    /* skipping first line because it contains field names */
    std::getline(in, line);

    /* read all stations */
    while (std::getline(in, line))
    {
        std::list<std::string> fields = BaseLib::splitString(line, '\t');

        if (fields.size() >= 5)
        {
            if (*fields.begin() == cName)  // add new layer
            {
                it = fields.begin();
                auto const pnt = parsePointCoordinates(it);

                // check if current layer has a thickness of 0.0.
                // if so skip it since it will mess with the vtk-visualisation
                // later on!
                if (pnt[2] != depth)
                {
                    if (newBorehole == nullptr)
                        OGS_FATAL("Trying to access a nullptr.");
                    newBorehole->addSoilLayer(pnt[0], pnt[1], pnt[2], sName);
                    sName = (*(++it));
                    depth = pnt[2];
                }
                else
                {
                    WARN(
                        "GMSInterface::readBoreholeFromGMS(): Skipped layer "
                        "'{:s}' in borehole '{:s}' because of thickness 0.0.",
                        sName, cName);
                }
            }
            else  // add new borehole
            {
                if (newBorehole != nullptr)
                {
                    newBorehole->setDepth((*newBorehole)[2] - depth);
                    boreholes.push_back(newBorehole);
                }
                cName = *fields.begin();
                it = fields.begin();
                auto const pnt = parsePointCoordinates(it);
                sName = (*(++it));
                newBorehole = GeoLib::StationBorehole::createStation(
                    cName, pnt[0], pnt[1], pnt[2], 0);
                depth = pnt[2];
            }
        }
        else
        {
            ERR("GMSInterface::readBoreholeFromGMS(): Error reading format.");
        }
    }
    // write the last borehole from the file
    if (newBorehole != nullptr)
    {
        newBorehole->setDepth((*newBorehole)[2] - depth);
        boreholes.push_back(newBorehole);
    }
    in.close();

    if (boreholes.empty())
    {
        return 0;
    }
    return 1;
}

void GMSInterface::writeBoreholesToGMS(
    const std::vector<GeoLib::Point*>* stations, const std::string& filename)
{
    std::ofstream out(filename.c_str(), std::ios::out);

    // write header
    out << "name"
        << "\t" << std::fixed << "X"
        << "\t"
        << "Y"
        << "\t"
        << "Z"
        << "\t"
        << "soilID"
        << "\n";

    for (auto station_as_point : *stations)
    {
        auto const* station =
            static_cast<GeoLib::StationBorehole*>(station_as_point);
        std::vector<GeoLib::Point*> const& profile = station->getProfile();
        std::vector<std::string> const& soilNames = station->getSoilNames();
        std::string current_soil_name;

        std::size_t nLayers = profile.size();
        for (std::size_t i = 1; i < nLayers; i++)
        {
            if ((i > 1) && (soilNames[i] == soilNames[i - 1]))
            {
                continue;
            }
            current_soil_name = soilNames[i];

            out << station->getName() << "\t" << std::fixed
                << (*(profile[i - 1]))[0] << "\t" << (*(profile[i - 1]))[1]
                << "\t" << (*(profile[i - 1]))[2] << "\t"
                << current_soil_name /*idx*/ << "\n";
        }
        out << station->getName() << "\t" << std::fixed
            << (*(profile[nLayers - 1]))[0] << "\t"
            << (*(profile[nLayers - 1]))[1] << "\t"
            << (*(profile[nLayers - 1]))[2] << "\t" << current_soil_name
            << "\n";  // this line marks the end of the borehole
    }

    out.close();
}

MeshLib::Mesh* GMSInterface::readMesh(const std::string& filename)
{
    std::string line;

    std::ifstream in(filename.c_str());
    if (!in.is_open())
    {
        ERR("GMSInterface::readMesh(): Could not open file {:s}.", filename);
        return nullptr;
    }

    // Read data from file
    std::getline(in, line);  // "MESH3D"
    if (line != "MESH3D" && line != "MESH2D")
    {
        ERR("GMSInterface::readMesh(): Could not read expected file "
            "header.");
        return nullptr;
    }
    bool const is_3d = (line == "MESH3D");

    std::string mesh_name = BaseLib::extractBaseNameWithoutExtension(filename);
    INFO("Reading SMS/GMS mesh...");
    std::vector<MeshLib::Node*> nodes;
    std::vector<MeshLib::Element*> elements;
    std::vector<int> mat_ids;
    std::map<unsigned, unsigned> id_map;

    // elements are listed before nodes in 3dm-format, therefore
    // traverse file twice and read first nodes and then elements
    std::string dummy;
    unsigned id(0);
    unsigned count(0);
    double x[3];
    // read nodes
    while (std::getline(in, line))
    {
        if (line[0] == 'N')  // "ND" for Node
        {
            std::stringstream str(line);
            str >> dummy >> id >> x[0] >> x[1] >> x[2];
            auto* node = new MeshLib::Node(x, id);
            id_map.insert(std::pair<unsigned, unsigned>(id, count++));
            nodes.push_back(node);
        }
    }
    in.close();

    // NOTE: Element types E8H (Hex), E4Q (Quad) are not implemented yet
    // read elements
    in.open(filename.c_str());
    std::getline(in, line);  // "MESH2D" / "MESH3D"
    unsigned node_idx[6];
    int mat_id(0);
    while (std::getline(in, line))
    {
        std::string element_id(line.substr(0, 3));
        std::stringstream str(line);

        if (element_id == "MES")  // "MESHNAME"
        {
            str >> dummy >> mesh_name;
            mesh_name = mesh_name.substr(1, mesh_name.length() - 2);
        }
        else if (!is_3d && element_id == "E3T")  // Triangle
        {
            str >> dummy >> id >> node_idx[0] >> node_idx[1] >> node_idx[2] >>
                mat_id;
            std::array<MeshLib::Node*, 3> tri_nodes;
            for (unsigned k(0); k < 3; k++)
            {
                tri_nodes[k] = nodes[id_map.find(node_idx[k])->second];
            }
            elements.push_back(new MeshLib::Tri(tri_nodes));
            mat_ids.push_back(mat_id);
        }
        else if (!is_3d && element_id == "E6T")  // Triangle
        {
            str >> dummy >> id >> node_idx[0] >> node_idx[3] >> node_idx[1] >>
                node_idx[4] >> node_idx[2] >> node_idx[5] >> mat_id;
            std::array<MeshLib::Node*, 3> tri_nodes;
            for (unsigned k(0); k < 3; k++)
            {
                tri_nodes[k] = nodes[id_map.find(node_idx[k])->second];
            }
            elements.push_back(new MeshLib::Tri(tri_nodes));
            mat_ids.push_back(mat_id);
        }
        else if (is_3d && element_id == "E6W")  // Prism
        {
            str >> dummy >> id >> node_idx[0] >> node_idx[1] >> node_idx[2] >>
                node_idx[3] >> node_idx[4] >> node_idx[5] >> mat_id;
            std::array<MeshLib::Node*, 6> prism_nodes;
            for (unsigned k(0); k < 6; k++)
            {
                prism_nodes[k] = nodes[id_map.find(node_idx[k])->second];
            }
            elements.push_back(new MeshLib::Prism(prism_nodes));
            mat_ids.push_back(mat_id);
        }
        else if (is_3d && element_id == "E4T")  // Tet
        {
            str >> dummy >> id >> node_idx[0] >> node_idx[1] >> node_idx[2] >>
                node_idx[3] >> mat_id;
            std::array<MeshLib::Node*, 4> tet_nodes;
            for (unsigned k(0); k < 4; k++)
            {
                tet_nodes[k] = nodes[id_map.find(node_idx[k])->second];
            }
            elements.push_back(new MeshLib::Tet(tet_nodes));
            mat_ids.push_back(mat_id);
        }
        // Pyramid (two versions exist for some reason)
        else if (is_3d && (element_id == "E4P" || element_id == "E5P"))
        {
            str >> dummy >> id >> node_idx[0] >> node_idx[1] >> node_idx[2] >>
                node_idx[3] >> node_idx[4] >> mat_id;
            std::array<MeshLib::Node*, 5> pyramid_nodes;
            for (unsigned k(0); k < 5; k++)
            {
                pyramid_nodes[k] = nodes[id_map.find(node_idx[k])->second];
            }
            elements.push_back(new MeshLib::Pyramid(pyramid_nodes));
            mat_ids.push_back(mat_id);
        }
        else if (element_id == "ND ")
        {              // Node
            continue;  // skip because nodes have already been read
        }
        else  // default
        {
            WARN(
                "GMSInterface::readMesh() - Element type '{:s}' not "
                "recognised.",
                element_id);
            return nullptr;
        }
    }

    in.close();
    INFO("finished.");

    MeshLib::Properties properties;
    if (mat_ids.size() == elements.size())
    {
        auto* const opt_pv = properties.createNewPropertyVector<int>(
            "MaterialIDs", MeshLib::MeshItemType::Cell);
        if (!opt_pv)
        {
            ERR("Could not create PropertyVector for material ids.");
            BaseLib::cleanupVectorElements(nodes, elements);
            return nullptr;
        }
        opt_pv->reserve(mat_ids.size());
        std::copy(mat_ids.cbegin(), mat_ids.cend(),
                  std::back_inserter(*opt_pv));
    }
    else
    {
        ERR("Ignoring Material IDs information (does not match number of "
            "elements).");
    }
    return new MeshLib::Mesh(mesh_name, nodes, elements,
                             true /* compute_element_neighbors */, properties);
}

}  // namespace FileIO
