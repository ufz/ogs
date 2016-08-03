/**
 * @copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "SwmmInterface.h"

#include <memory>
#include <utility>

#include "BaseLib/FileTools.h"
#include "GeoLib/GEOObjects.h"
#include "GeoLib/Point.h"
#include "GeoLib/Polyline.h"
#include "GeoLib/Polygon.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/Elements/Line.h"
#include "MeshLib/Properties.h"

#include "Applications/ApplicationsLib/LogogSetup.h"


SwmmInterface::SwmmInterface(std::string const& swmm_base_name)
: _base_name (swmm_base_name)
{
}

SwmmInterface::SwmmInterface(MeshLib::Mesh* mesh, std::string const& swmm_base_name)
: _base_name (swmm_base_name), _mesh (mesh)
{
}

int SwmmInterface::isSwmmInputFile(std::string const& inp_file_name)
{
    std::ifstream in ( inp_file_name.c_str() );
    int ret = isSwmmFile(inp_file_name, in, "inp");

    if (ret != 0)
        return ret;

    std::string line;
    bool header_found (false);
    std::size_t pos_beg (0), pos_end (0);
    while (!header_found)
    {
        getline(in, line);
        pos_beg = line.find_first_not_of(' ', pos_end);
        pos_end = line.find_first_of(" \n", pos_beg);

        if (line.empty() || pos_beg==pos_end || line.compare(pos_beg,2,";;") == 0)
            continue;

        if (line == "[TITLE]")
            header_found = true;
        else
        {
            ERR ("SWMMInterface: input file type %s not recognised.",
                BaseLib::getFileExtension(inp_file_name).c_str());
            return -3;
        }
    }

    in.close();
    return ret;
}

int SwmmInterface::isSwmmOutputFile(std::string const& out_file_name)
{
    std::ifstream in ( out_file_name.c_str() );
    int ret = isSwmmFile(out_file_name, in, "out");
    in.close();
    return ret;
}

int SwmmInterface::isSwmmFile(std::string const& inp_file_name, std::ifstream &in, std::string const& ext)
{
    if (BaseLib::getFileExtension(inp_file_name) != ext)
    {
        ERR ("SWMMInterface: file extension %s not recognised (should be *.%s).",
            BaseLib::getFileExtension(inp_file_name).c_str(), ext.c_str());
        return -2;
    }

    if (!in.is_open())
    {
        ERR ("SWMMInterface: Could not open input file %s.", inp_file_name.c_str());
        return -1;
    }
    return 0;
}

int SwmmInterface::readPolygons(std::ifstream &in, std::vector<GeoLib::Polyline*> &lines,
    std::vector<std::string> &line_names, std::vector<GeoLib::Point*> &points,
    std::vector<std::string> &pnt_names)
{
    bool finished (false);
    std::size_t id (points.size());
    std::string line;
    std::string line_name("");
    GeoLib::Polyline* p (nullptr);
    while (!finished)
    {
        getline(in, line);
        std::size_t pos_end = 0;
        std::size_t pos_beg = line.find_first_not_of(' ', pos_end);
        pos_end = line.find_first_of(" \n", pos_beg);

        if (line.empty() || pos_beg==pos_end)
        {
            if (p != nullptr)
                lines.push_back(p);
            finished = true;
            return 0;
        }

        if (line.compare(pos_beg,2,";;") == 0)
            continue;

        pos_beg = line.find_first_not_of(' ', 0);
        pos_end = line.find_first_of(' ', pos_beg);
        if (line.substr(pos_beg, pos_end-pos_beg) != line_name)
        {
            if (p != nullptr)
                lines.push_back(p);

            line_name = line.substr(pos_beg, pos_end-pos_beg);
            p = new GeoLib::Polyline(points);
            line_names.push_back(line_name);
        }

        GeoLib::Point* pnt = new GeoLib::Point(0, 0, 0, id);
        pnt_names.push_back("");
        for (std::size_t i=0; i<2; ++i)
        {
            pos_beg = line.find_first_not_of(' ', pos_end);
            pos_end = line.find_first_of(' ', pos_beg);
            if (pos_beg != std::string::npos)
                (*pnt)[i] = BaseLib::str2number<double>(line.substr(pos_beg, pos_end-pos_beg));
            else
            {
                ERR("SwmmInterface::readPolygons(): error reading coordinate %d of object %s.", 
                    i, line_names[line_names.size()-1].c_str());
                return false;
            }
        }
        points.push_back(pnt);
        id++;
        p->addPoint(points.size()-1);
    }
    return 0;
}

int SwmmInterface::SwmmInputToGeometry(std::string const& inp_file_name, GeoLib::GEOObjects &geo_objects)
{
    std::ifstream in ( inp_file_name.c_str() );
    int const ret (isSwmmFile(inp_file_name, in, "inp"));
    if ( ret != 0)
        return ret;

    std::unique_ptr<std::vector<GeoLib::Point*>> points (new std::vector<GeoLib::Point*>);
    std::unique_ptr<std::vector<GeoLib::Polyline*>> lines (new std::vector<GeoLib::Polyline*>);
    std::vector<std::string> pnt_names;
    std::vector<std::string> line_names;

    std::string geo_name = BaseLib::extractBaseNameWithoutExtension(inp_file_name);
    std::string line;
    while ( getline(in, line) )
    {
        if (line == "[COORDINATES]")
        {
            if (readCoordinates<GeoLib::Point>(in, *points, pnt_names) != 0)
                return -1;
        }
        if (line == "[VERTICES]")
        {
            if (readCoordinates<GeoLib::Point>(in, *points, pnt_names) != 0)
                return -1;
        }
        if (line == "[Polygons]")
        {
            if (readPolygons(in, *lines, line_names, *points, pnt_names) != 0)
                return -2;
        }
        if (line == "[SYMBOLS]")
        {
            if (readCoordinates<GeoLib::Point>(in, *points, pnt_names) != 0)
                return -1;
        }
    }

    if (!points->empty())
    {
        if (points->size() != pnt_names.size())
        {
            ERR ("Lengt of point vector and point name vector do not match.");
            return -3;
        }
        std::map<std::string, std::size_t> *name_id_map (new std::map<std::string, std::size_t>);
        std::size_t const n_names (pnt_names.size());
        for (std::size_t i=0; i<n_names; ++i)
        {
            if (pnt_names[i] != "")
                name_id_map->insert(std::make_pair(pnt_names[i], i));
        }
        geo_objects.addPointVec(std::move(points), geo_name, name_id_map);
    }
    else
    {
        ERR ("No points found in file");
        return -4;
    }
    if (!lines->empty())
    {
        if (lines->size() != line_names.size())
        {
            ERR ("Lengt of line vector and line name vector do not match.");
            return -5;
        }
        std::map<std::string, std::size_t> *name_id_map (new std::map<std::string, std::size_t>);
        std::size_t const n_names (line_names.size());
        for (std::size_t i=0; i<n_names; ++i)
            name_id_map->insert(std::make_pair(line_names[i], i));
        geo_objects.addPolylineVec(std::move(lines), geo_name, name_id_map);
    }

    return 0;
}

int SwmmInterface::readNodeData(std::ifstream &in, std::vector<MeshLib::Node*> &nodes, std::map<std::string,
    std::size_t> const& name_id_map, std::vector<double> &max_depth, bool read_max_depth)
{
    bool finished (false);
    std::string line;
    while (!finished)
    {
        getline(in, line);
        std::size_t pos_end = 0;
        std::size_t pos_beg = line.find_first_not_of(' ', pos_end);
        pos_end = line.find_first_of(" \n", pos_beg);

        if (line.empty() || pos_beg==pos_end)
        {
            finished = true;
            return 0;
        }

        if (line.compare(pos_beg,2,";;") == 0)
            continue;

        pos_beg = line.find_first_not_of(' ', 0);
        pos_end = line.find_first_of(' ', pos_beg);
        std::string const current_name (line.substr(pos_beg, pos_end-pos_beg));
        auto it = name_id_map.find(current_name);
        if (it == name_id_map.end())
        {
            ERR ("SwmmInterface::readNodeData(): Name %s not found coordinates map.", current_name.c_str());
            return -1;
        }
        std::size_t id = it->second;
        pos_beg = line.find_first_not_of(' ', pos_end);
        pos_end = line.find_first_of(' ', pos_beg);
        if (pos_beg != std::string::npos)
            (*nodes[id])[2] = BaseLib::str2number<double>(line.substr(pos_beg, pos_end-pos_beg));
        else
        {
            ERR ("SwmmInterface::readNodeData(): error reading elevation of object %s.", current_name.c_str());
            return -2;
        }

        if (!read_max_depth)
        {
            max_depth[id] = 0;
            continue;
        }

        pos_beg = line.find_first_not_of(' ', pos_end);
        pos_end = line.find_first_of(' ', pos_beg);
        if (pos_beg != std::string::npos)
            max_depth[id] = BaseLib::str2number<double>(line.substr(pos_beg, pos_end-pos_beg));
        else
        {
            ERR ("SwmmInterface::readNodeData(): error reading max. depth of object %s.", current_name.c_str());
            return -3;
        }
    }
    return 0;
}

int SwmmInterface::readLineElements(std::ifstream &in, std::vector<MeshLib::Element*> &elements,
    std::vector<MeshLib::Node*> const& nodes, std::map<std::string, std::size_t> const& name_id_map)
{
    bool finished (false);
    std::string line;
    while (!finished)
    {
        getline(in, line);
        std::size_t pos_end = 0;
        std::size_t pos_beg = line.find_first_not_of(' ', pos_end);
        pos_end = line.find_first_of(" \n", pos_beg);

        if (line.empty() || pos_beg==pos_end)
        {
            finished = true;
            return 0;
        }

        if (line.compare(pos_beg,2,";;") == 0)
            continue;

        pos_beg = line.find_first_not_of(' ', 0);
        pos_end = line.find_first_of(' ', pos_beg);
        // the name of the conduit is currently not needed but might be relevant later
        std::string const current_name (line.substr(pos_beg, pos_end-pos_beg));

        pos_beg = line.find_first_not_of(' ', pos_end);
        pos_end = line.find_first_of(' ', pos_beg);
        std::string const inlet (line.substr(pos_beg, pos_end-pos_beg));
        auto i_it = name_id_map.find(inlet);
        if (i_it == name_id_map.end())
        {
            ERR ("SwmmInterface::readLineElements(): Inlet node %s not found coordinates map.", inlet.c_str());
            return -1;
        }

        pos_beg = line.find_first_not_of(' ', pos_end);
        pos_end = line.find_first_of(' ', pos_beg);
        std::string const outlet (line.substr(pos_beg, pos_end-pos_beg));
        auto o_it = name_id_map.find(outlet);
        if (o_it == name_id_map.end())
        {
            ERR ("SwmmInterface::readLineElements(): Outlet node %s not found coordinates map.", outlet.c_str());
            return -2;
        }

        //MeshLib::Node lines_nodes[2]
        std::array<MeshLib::Node*, 2> line_nodes = { nodes[i_it->second], nodes[o_it->second] };
        elements.push_back(new MeshLib::Line(line_nodes));
    }
    return 0;
}

MeshLib::Mesh* SwmmInterface::SwmmInputToLineMesh()
{
    std::string const inp_file_name (_base_name + ".inp");
    std::ifstream in ( inp_file_name.c_str() );
    int const ret (isSwmmFile(inp_file_name, in, "inp"));
    if ( ret != 0)
        return nullptr;

    std::vector< MeshLib::Node* > nodes;
    std::vector<std::string> id_name_map;
    std::string line;
    while ( getline(in, line) )
    {
        if (line == "[COORDINATES]")
        {
            INFO ("Reading coordinates...");
            //INFO ("Reading coordinates...");
            if (readCoordinates<MeshLib::Node>(in, nodes, id_name_map) != 0)
                return nullptr;
        }
        if (line == "[VERTICES]")
        {
            INFO ("Reading vertices...");
            //INFO ("Reading vertices...");
            if (readCoordinates(in, nodes, id_name_map) != 0)
                return nullptr;
        }
    }

    // After end of file is reached, create name-id-map and
    // start reading again to get line elements and node data.
    std::map< std::string, std::size_t> name_id_map;
    std::size_t n_nodes (nodes.size());
    for (std::size_t i=0; i<n_nodes; ++i)
        name_id_map[id_name_map[i]] = i;
    in.clear();
    in.seekg(0, in.beg);

    std::vector< MeshLib::Element* > elements;
    std::vector<double> max_depth;
    max_depth.resize(n_nodes);
    std::size_t const n_types = 3;
    std::array< std::size_t, n_types> n_elem_types;
    while ( getline(in, line) )
    {
        if (line == "[JUNCTIONS]")
        {
            INFO ("Reading junctions...");
            if (readNodeData(in, nodes, name_id_map, max_depth, true) != 0)
                return nullptr;
        }
        if (line == "[OUTFALLS]")
        {
            INFO ("Reading outfalls...");
            if (readNodeData(in, nodes, name_id_map, max_depth, false) != 0)
                return nullptr;
        }
        if (line == "[STORAGE]")
        {
            INFO ("Reading storages...");
            if (readNodeData(in, nodes, name_id_map, max_depth, true) != 0)
                return nullptr;
        }
        if (line == "[CONDUITS]")
        {
            INFO ("Reading conduits...");
            if (readLineElements(in, elements, nodes, name_id_map) != 0)
                return nullptr;
            n_elem_types[0] = elements.size();
        }
        if (line == "[PUMPS]")
        {
            INFO ("Reading pumps...");
            if (readLineElements(in, elements, nodes, name_id_map) != 0)
                return nullptr;
            n_elem_types[1] = elements.size();
        }
        if (line == "[WEIRS]")
        {
            INFO ("Reading weirs...")
            if (readLineElements(in, elements, nodes, name_id_map) != 0)
                return nullptr;
            n_elem_types[2] = elements.size();
        }
    }

    if (nodes.empty() || elements.empty())
        return nullptr;

    MeshLib::Properties props;
    boost::optional< MeshLib::PropertyVector<int>& > mat_ids =
        props.createNewPropertyVector<int>("MaterialIDs", MeshLib::MeshItemType::Cell, 1);
    mat_ids->resize(elements.size(), 0);
    for (std::size_t i=1; i<n_types; ++i)
        std::fill(mat_ids->begin()+n_elem_types[i-1], mat_ids->begin()+n_elem_types[i], i);

    if (nodes.size() == max_depth.size())
    {
        boost::optional< MeshLib::PropertyVector<double>& > depth =
            props.createNewPropertyVector<double>("Max Depth", MeshLib::MeshItemType::Node, 1);
        depth->reserve(max_depth.size());
        std::copy(max_depth.cbegin(), max_depth.cend(), std::back_inserter(*depth));
    }
    else
        ERR ("Size of max depth array does not fit number of elements. Skipping array.");

    return new MeshLib::Mesh(_base_name, nodes, elements, props);
}
