/**
 * @copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "SwmmInterface.h"

#include <cctype>
#include <iostream>
#include <utility>

#include "BaseLib/FileTools.h"
#include "BaseLib/StringTools.h"
#include "GeoLib/GEOObjects.h"
#include "GeoLib/Point.h"
#include "GeoLib/Polyline.h"
#include "GeoLib/Polygon.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/Elements/Line.h"
#include "MeshLib/Properties.h"

#include "Applications/FileIO/CsvInterface.h"

#include "Applications/ApplicationsLib/LogogSetup.h"

#include "ThirdParty/SWMMInterface/swmm5_iface.h"


namespace FileIO
{

/// Variables that always exist for subcatchments. There might be more.
const std::array<std::string,9> subcatchment_vars =
{
    "rainfall",
    "snow depth",
    "evaporation loss",
    "infiltration losses",
    "runoff rate",
    "groundwater outflow",
    "groundwater head",
    "moisture content",
    "concentration of pollutant"
};

/// Variables that always exist for nodes. There might be more.
const std::array<std::string,7> node_vars =
{
    "water depth",
    "hydraulic head",
    "volume of stored water",
    "lateral inflow",
    "total inflow",
    "flow lost to flooding",
    "concentration of pollutant"
};

/// Variables that always exist for links. There might be more.
const std::array<std::string,6> link_vars =
{
    "flow rate",
    "flow depth",
    "flow velocity",
    "flow volume",
    "fraction conduit/non-conduit",
    "concentration of pollutant"
};

/// All variables that exist for the system.
const std::array<std::string,15> system_vars =
{
    "air temperature",
    "rainfall",
    "snow depth",
    "evaporation + infiltration loss",
    "runoff flow",
    "dry weather inflow",
    "groundwater inflow",
    "RDII inflow",
    "direct inflow",
    "total lateral inflow",
    "flow lost to flooding",
    "flow leaving through outfalls",
    "volume of stored water",
    "actual evaporation rate",
    "potential evaporation rate"
};

/// Number of base variables for the four object types (subcatchments/nodes/links/system).
std::array<std::size_t,4> const n_obj_params = { 8, 6, 5, 15 };

std::unique_ptr<SwmmInterface> SwmmInterface::create(std::string const& file_name)
{
    //The input/output methods take a base name and check if the corresponding i/o file for that base name exists.
    //This check takes any swmm project file, i.e. [base name].[extension] which needs to be at least 5 chars
    //because the extension is always 3 chars.
    if (file_name.length() < 5)
        return nullptr;

    if (!SwmmInterface::isSwmmInputFile(file_name))
        return nullptr;

    std::string const base_name (file_name.substr(0, file_name.length() - 4));
    SwmmInterface* swmm = new SwmmInterface(base_name);
    if (swmm->readSwmmInputToLineMesh())
        return std::unique_ptr<SwmmInterface>(swmm);

    ERR ("Error creating mesh from SWMM file.");
    delete swmm;
    return nullptr;
}

SwmmInterface::SwmmInterface(std::string const& swmm_base_name)
: _base_name (swmm_base_name), _mesh(nullptr)
{
}

SwmmInterface::~SwmmInterface()
{
    for (Subcatchment& sc  : _subcatchments)
        delete sc.outline;

    for (GeoLib::Point* pnt : _subcatchment_points)
        delete pnt;
}

bool SwmmInterface::isSwmmInputFile(std::string const& inp_file_name)
{
    std::string path (inp_file_name);
    std::transform(path.begin(), path.end(), path.begin(), std::tolower);
    if (BaseLib::getFileExtension(path) != "inp")
    {
        ERR ("SWMMInterface: %s is not a SWMM input file.", inp_file_name.c_str());
        return false;
    }

    std::ifstream in ( inp_file_name.c_str() );
    if (!in.is_open())
    {
        ERR ("SWMMInterface: Could not open input file %s.", inp_file_name.c_str());
        return false;
    }

    std::string line;
    bool header_found (false);
    std::size_t pos_beg (0), pos_end (0);
    while (!header_found)
    {
        if (!std::getline(in, line))
            return false;

        pos_beg = line.find_first_not_of(' ', pos_end);
        pos_end = line.find_first_of(" \n", pos_beg);

        // skip empty or comment lines at the beginning of the file
        if (line.empty() || pos_beg==pos_end || isCommentLine(line))
            continue;

        if (line == "[TITLE]")
            header_found = true;
        else
        {
            INFO ("SWMMInterface: input file type %s not recognised.",
                BaseLib::getFileExtension(inp_file_name).c_str());
            return false;
        }
    }

    in.close();
    return true;
}

bool SwmmInterface::existsSwmmOutputFile() const
{
    std::string const outfile (_base_name + ".out");
    if (OpenSwmmOutFile(const_cast<char*>(outfile.c_str())) != 0)
        return false;
    return true;
}

template <typename T>
bool SwmmInterface::readCoordinates(std::ifstream &in, std::vector<T*> &points, std::vector<std::string> &names)
{
    std::size_t id (points.size());
    std::string line;

    while (std::getline(in, line))
    {
        if (isSectionFinished(line))
            return true;

        if (isCommentLine(line))
            continue;

        std::vector<std::string> split_str (BaseLib::splitString(line));
        if (split_str.size() != 3)
        {
            ERR ("Format not recognised.");
            return false;
        }
        names.push_back(split_str[0]);

        double const x = BaseLib::str2number<double>(split_str[1]);
        double const y = BaseLib::str2number<double>(split_str[2]);
        T* pnt = new T(x, y, 0, id);
        points.push_back(pnt);
        id++;
    }
    return true;
}

bool SwmmInterface::readPolygons(std::ifstream &in, std::vector<GeoLib::Polyline*> &lines,
    std::vector<std::string> &ply_names, std::vector<GeoLib::Point*> &points,
    std::vector<std::string> &pnt_names)
{
    bool finished (false);
    std::size_t id (points.size());
    std::string line;
    std::string polygon_name("");
    GeoLib::Polyline* p (nullptr);
    while (std::getline(in, line))
    {
        if (isSectionFinished(line))
            break;

        if (isCommentLine(line))
            continue;

        std::vector<std::string> split_str (BaseLib::splitString(line));
        if (split_str.size() != 3)
        {
            ERR ("Polygon format not recognised.");
            delete p;
            return false;
        }

        // if a new polygon starts, add the old one to the vector
        if (split_str[0] != polygon_name)
        {
            if (p != nullptr)
                lines.push_back(p);

            polygon_name = split_str[0];
            p = new GeoLib::Polyline(points);
            ply_names.push_back(polygon_name);
        }

        double const x = BaseLib::str2number<double>(split_str[1]);
        double const y = BaseLib::str2number<double>(split_str[2]);
        points.push_back(new GeoLib::Point(x, y, 0, id));
        p->addPoint(points.size()-1);
        pnt_names.push_back("");
        id++;
    }

    // when the section is finished, add the last polygon
    if (p != nullptr)
        lines.push_back(p);

    return true;
}

bool SwmmInterface::addPointElevation(std::ifstream &in,
    std::vector<GeoLib::Point*> &points, std::map<std::string,
    std::size_t> const& name_id_map)
{
    std::string line;
    while (std::getline(in, line))
    {
        if (isSectionFinished(line))
            return true;

        if (isCommentLine(line))
            continue;

        std::vector<std::string> const split_str (BaseLib::splitString(line));
        // Junctions = 6, Outfalls = 4, Storage = 8
        if (split_str.size() < 4)
        {
            ERR ("Format not recognised.");
            return false;
        }
        std::string const current_name (split_str[0]);
        auto const it = name_id_map.find(current_name);
        if (it == name_id_map.end())
        {
            ERR ("SwmmInterface::addPointElevation(): Name %s not found in coordinates map.", current_name.c_str());
            return false;
        }
        std::size_t const id = it->second;
        (*points[id])[2] = BaseLib::str2number<double>(split_str[1]);
    }
    return true;
}

bool SwmmInterface::readLinksAsPolylines(std::ifstream &in,
    std::vector<GeoLib::Polyline*> &lines, std::vector<std::string> &line_names,
    std::vector<GeoLib::Point*> const& points, std::map<std::string, std::size_t> const& point_names)
{
    std::string line;
    while (std::getline(in, line))
    {
        if (isSectionFinished(line))
            return true;

        if (isCommentLine(line))
            continue;

        std::vector<std::string> const split_str (BaseLib::splitString(line));
        // Conduits = 9, Pumps = 7, Weirs = 8
        if (split_str.size() < 7)
        {
            ERR ("Conduit format not recognised.");
            return false;
        }

        std::string const inlet (split_str[1]);
        auto const i_it = point_names.find(inlet);
        if (i_it == point_names.end())
        {
            ERR ("SwmmInterface::readLineElements(): Inlet node %s not found in coordinates map.", inlet.c_str());
            return false;
        }

        std::string const outlet (split_str[2]);
        auto const o_it = point_names.find(outlet);
        if (o_it == point_names.end())
        {
            ERR ("SwmmInterface::readLineElements(): Outlet node %s not found in coordinates map.", outlet.c_str());
            return false;
        }
        GeoLib::Polyline* ply = new GeoLib::Polyline(points);
        std::size_t a (i_it->second);
        ply->addPoint(i_it->second);
        ply->addPoint(o_it->second);
        lines.push_back(ply);
        line_names.push_back(split_str[0]);
    }
    return true;
}

/// Deletes the geometric objects and returns false
bool geometryCleanup(std::vector<GeoLib::Point*> &points, std::vector<GeoLib::Polyline*> &lines)
{
    for (auto line : lines)
        delete line;
    for (auto point : points)
        delete point;
    return false;
}

bool SwmmInterface::convertSwmmInputToGeometry(std::string const& inp_file_name,
    GeoLib::GEOObjects &geo_objects, bool add_subcatchments)
{
    if (!isSwmmInputFile(inp_file_name))
        return false;

    std::ifstream in ( inp_file_name.c_str() );
    if (!in.is_open())
    {
        ERR ("SWMMInterface: Could not open input file %s.", inp_file_name.c_str());
        return false;
    }

    std::unique_ptr<std::vector<GeoLib::Point*>> points (new std::vector<GeoLib::Point*>);
    std::unique_ptr<std::vector<GeoLib::Polyline*>> lines (new std::vector<GeoLib::Polyline*>);
    std::vector<std::string> pnt_names;
    std::vector<std::string> line_names;

    std::string geo_name = BaseLib::extractBaseNameWithoutExtension(inp_file_name);
    std::string line;
    while ( std::getline(in, line) )
    {
        if (line == "[COORDINATES]")
        {
            if (!readCoordinates<GeoLib::Point>(in, *points, pnt_names))
                return geometryCleanup(*points, *lines);
        }
        if (line == "[VERTICES]")
        {
            if (!readCoordinates<GeoLib::Point>(in, *points, pnt_names))
                return geometryCleanup(*points, *lines);
        }
        if (line == "[Polygons]" && add_subcatchments)
        {
            if (!readPolygons(in, *lines, line_names, *points, pnt_names))
                return geometryCleanup(*points, *lines);
        }
        if (line == "[SYMBOLS]")
        {
            if (!readCoordinates<GeoLib::Point>(in, *points, pnt_names))
                return geometryCleanup(*points, *lines);
        }
    }

    if (points->empty())
    {
        ERR ("No points found in file");
        return false;
    }
    if (points->size() != pnt_names.size())
    {
        ERR ("Length of point vector and point name vector do not match.");
        return geometryCleanup(*points, *lines);
    }

    std::map<std::string, std::size_t> *name_id_map (new std::map<std::string, std::size_t>);
    std::size_t const n_names (pnt_names.size());
    for (std::size_t i=0; i<n_names; ++i)
    {
        if (!pnt_names[i].empty())
            name_id_map->insert(std::make_pair(pnt_names[i], i));
    }

    // rewind stream and read links between junctions
    in.clear();
    in.seekg(0, in.beg);

    while ( std::getline(in, line) )
    {
        if (line == "[JUNCTIONS]")
        {
            INFO ("Reading point elevation...");
            if (!addPointElevation(in, *points, *name_id_map))
                return geometryCleanup(*points, *lines);
        }
        if (line == "[CONDUITS]")
        {
            INFO ("Reading conduits...");
            if (!readLinksAsPolylines(in, *lines, line_names, *points, *name_id_map))
                return geometryCleanup(*points, *lines);
        }
        else if (line == "[PUMPS]")
        {
            INFO ("Reading pumps...");
            if (!readLinksAsPolylines(in, *lines, line_names, *points, *name_id_map))
                return geometryCleanup(*points, *lines);
        }
        else if (line == "[WEIRS]")
        {
            INFO ("Reading weirs...")
            if (!readLinksAsPolylines(in, *lines, line_names, *points, *name_id_map))
                return geometryCleanup(*points, *lines);
        }
    }

    geo_objects.addPointVec(std::move(points), geo_name, name_id_map);
    if (!lines->empty())
    {
        if (lines->size() != line_names.size())
        {
            ERR ("Lengt of line vector and line name vector do not match.");
            geo_objects.removePointVec(geo_name);
            for (auto ply : *lines)
                delete ply;
            return false;
        }
        std::map<std::string, std::size_t> *line_id_map (new std::map<std::string, std::size_t>);
        std::size_t const n_names (line_names.size());
        for (std::size_t i=0; i<n_names; ++i)
            line_id_map->insert(std::make_pair(line_names[i], i));
        geo_objects.addPolylineVec(std::move(lines), geo_name, line_id_map);
    }
    return true;
}

bool SwmmInterface::readNodeData(std::ifstream &in, std::vector<MeshLib::Node*> &nodes, std::map<std::string,
    std::size_t> const& name_id_map, std::vector<double> &max_depth, bool read_max_depth)
{
    std::string line;
    while (std::getline(in, line))
    {
        if (isSectionFinished(line))
            return true;

        if (isCommentLine(line))
            continue;

        std::vector<std::string> const split_str (BaseLib::splitString(line));
        // Junctions = 6, Outfalls = 4, Storage = 8
        if (split_str.size() < 3)
        {
            ERR ("Format not recognised.");
            return false;
        }
        std::string const current_name (split_str[0]);
        auto const it = name_id_map.find(current_name);
        if (it == name_id_map.end())
        {
            ERR ("SwmmInterface::readNodeData(): Name %s not found in coordinates map.", current_name.c_str());
            return false;
        }
        std::size_t const id = it->second;
        (*nodes[id])[2] = BaseLib::str2number<double>(split_str[1]);

        if (read_max_depth)
            max_depth[id] = BaseLib::str2number<double>(split_str[2]);
        else
            max_depth[id] = 0;
    }
    return true;
}

bool SwmmInterface::readLineElements(std::ifstream &in, std::vector<MeshLib::Element*> &elements,
    std::vector<MeshLib::Node*> const& nodes, std::map<std::string, std::size_t> const& name_id_map)
{
    std::string line;
    while (std::getline(in, line))
    {
        if (isSectionFinished(line))
            return true;

        if (isCommentLine(line))
            continue;

        std::vector<std::string> const split_str (BaseLib::splitString(line));
        // Conduits = 9, Pumps = 7, Weirs = 8
        if (split_str.size() < 7)
        {
            ERR ("Conduit format not recognised.");
            return false;
        }

        std::string const inlet (split_str[1]);
        auto const i_it = name_id_map.find(inlet);
        if (i_it == name_id_map.end())
        {
            ERR ("SwmmInterface::readLineElements(): Inlet node %s not found in coordinates map.", inlet.c_str());
            return false;
        }

        std::string const outlet (split_str[2]);
        auto const o_it = name_id_map.find(outlet);
        if (o_it == name_id_map.end())
        {
            ERR ("SwmmInterface::readLineElements(): Outlet node %s not found in coordinates map.", outlet.c_str());
            return false;
        }

        std::array<MeshLib::Node*, 2> const line_nodes = { nodes[i_it->second], nodes[o_it->second] };
        elements.push_back(new MeshLib::Line(line_nodes));
        _id_linkname_map.push_back(split_str[0]);
    }
    return true;
}

bool SwmmInterface::readSubcatchments(std::ifstream &in, std::map< std::string, std::size_t> const& name_id_map)
{
    std::string line;
    while (getline(in, line))
    {
        if (isSectionFinished(line))
            return true;

        if (isCommentLine(line))
            continue;

        std::vector<std::string> const split_str (BaseLib::splitString(line));
        if (split_str.size() < 8)
        {
            ERR ("Subcatchment format not recognised.");
            return false;
        }

        Subcatchment sc;
        sc.name = split_str[0];
        sc.rain_gauge = std::numeric_limits<std::size_t>::max();
        std::size_t const n_gauges (_rain_gauges.size());
        for (std::size_t i=0; i<n_gauges; ++i)
        {
            if (_rain_gauges[i].first.getName() == split_str[1])
            {
                sc.rain_gauge = i;
                break;
            }
        }
        if (sc.rain_gauge == std::numeric_limits<std::size_t>::max())
        {
            ERR ("Rain gauge for subcatchment \"%s\" not found.", split_str[0].c_str());
            return false;
        }

        sc.outlet =  std::numeric_limits<std::size_t>::max();
        auto const it = name_id_map.find(split_str[2]);
        if (it == name_id_map.end())
        {
            ERR ("Outlet node for subcatchment \"%s\" not found.", split_str[0].c_str());
            return false;
        }
        sc.outlet = it->second;
        sc.area = BaseLib::str2number<double>(split_str[3]);
        _subcatchments.push_back(sc);
    }

    return true;
}

bool SwmmInterface::readSwmmInputToLineMesh()
{
    if (_mesh != nullptr)
    {
        ERR ("Mesh already exists.");
        return false;
    }

    std::string const inp_file_name (_base_name + ".inp");
    if (!isSwmmInputFile(inp_file_name))
        return false;

    std::ifstream in ( inp_file_name.c_str() );
    if (!in.is_open())
    {
        ERR ("SWMMInterface: Could not open input file %s.", inp_file_name.c_str());
        return false;
    }

    _id_nodename_map.clear();
    std::vector< MeshLib::Node* > nodes;
    std::string line;
    while ( getline(in, line) )
    {
        if (line == "[COORDINATES]")
        {
            INFO ("Reading coordinates...");
            if (!readCoordinates<MeshLib::Node>(in, nodes, _id_nodename_map))
                return false;
        }
        /* TODO: check if needed
        else if (line == "[VERTICES]")
        {
            INFO ("Reading vertices...");
            if (!readCoordinates(in, nodes, _id_nodename_map))
                return false;
        }
        */
        else if (line == "[SYMBOLS]")
        {
            INFO ("Reading symbols...");
            std::vector<GeoLib::Point*> points;
            std::vector<std::string> names;
            if (!readCoordinates(in,points, names))
                return false;
            for (std::size_t i=0; i<points.size(); ++i)
            {
                GeoLib::Station stn (points[i], names[i]);
                _rain_gauges.push_back(std::pair<GeoLib::Station, std::string>(stn, ""));
            }
        }
    }

    if (nodes.empty())
        return false;

    // After end of file is reached, create name-id-map and
    // start reading again to get line elements and node data.
    std::map< std::string, std::size_t> name_id_map;
    std::size_t const n_nodes (nodes.size());
    for (std::size_t i=0; i<n_nodes; ++i)
        name_id_map[_id_nodename_map[i]] = i;
    in.clear();
    in.seekg(0, in.beg);

    std::vector< MeshLib::Element* > elements;
    std::vector<double> max_depth(n_nodes);
    std::size_t const n_types = 3;
    std::array< std::size_t, n_types> n_elem_types {{0,0,0}};
    while ( getline(in, line) )
    {
        if (line == "[RAINGAGES]")
        {
            if (!_rain_gauges.empty())
                addRainGaugeTimeSeriesLocations(in);
        }
        else if (line == "[SUBCATCHMENTS]")
        {
            INFO ("Reading subcatchment information...");
            if (!readSubcatchments(in, name_id_map))
                return false;
        }
        else if (line == "[SUBAREAS]")
        {
            // more subcatchment variables, not yet implemented
        }
        else if (line == "[INFILTRATION]")
        {
            // more subcatchment variables, not yet implemented
        }
        else if (line == "[JUNCTIONS]")
        {
            INFO ("Reading junctions...");
            if (!readNodeData(in, nodes, name_id_map, max_depth, true))
                return false;
        }
        else if (line == "[OUTFALLS]")
        {
            INFO ("Reading outfalls...");
            if (!readNodeData(in, nodes, name_id_map, max_depth, false))
                return false;
        }
        else if (line == "[STORAGE]")
        {
            INFO ("Reading storages...");
            if (!readNodeData(in, nodes, name_id_map, max_depth, true))
                return false;
        }
        else if (line == "[CONDUITS]")
        {
            INFO ("Reading conduits...");
            if (!readLineElements(in, elements, nodes, name_id_map))
                return false;
            n_elem_types[0] = elements.size();
        }
        else if (line == "[PUMPS]")
        {
            INFO ("Reading pumps...");
            if (!readLineElements(in, elements, nodes, name_id_map))
                return false;
            n_elem_types[1] = elements.size();
        }
        else if (line == "[WEIRS]")
        {
            INFO ("Reading weirs...")
            if (!readLineElements(in, elements, nodes, name_id_map))
                return false;
            n_elem_types[2] = elements.size();
        }
        else if (line == "[POLLUTANTS]")
        {
            if (!readPollutants(in))
                return false;
        }
        else if (line == "[Polygons]")
        {
            INFO ("Reading subcatchments...");
            std::vector<GeoLib::Polyline*> lines;
            std::vector<std::string> line_names;
            std::vector<std::string> tmp_names; // polygon points are nameless but the method requires a vector
            if (!readPolygons(in, lines, line_names, _subcatchment_points, tmp_names))
                return false;

            if (!matchSubcatchmentsWithPolygons(lines, line_names))
                return false;
        }
    }

    if (elements.empty())
    {
        for (MeshLib::Node* node : nodes)
            delete node;
        return false;
    }

    MeshLib::Properties props;
    boost::optional< MeshLib::PropertyVector<int>& > mat_ids =
        props.createNewPropertyVector<int>("MaterialIDs", MeshLib::MeshItemType::Cell, 1);
    mat_ids->resize(elements.size(), 0);
    for (std::size_t i=1; i<n_types; ++i)
    {
        if (n_elem_types[i] > 0)
            std::fill(mat_ids->begin()+n_elem_types[i-1], mat_ids->begin()+n_elem_types[i], i);
    }

    if (nodes.size() == max_depth.size())
    {
        boost::optional< MeshLib::PropertyVector<double>& > depth =
            props.createNewPropertyVector<double>("Max Depth", MeshLib::MeshItemType::Node, 1);
        depth->reserve(max_depth.size());
        std::copy(max_depth.cbegin(), max_depth.cend(), std::back_inserter(*depth));
    }
    else
        ERR ("Size of max depth array does not fit number of elements. Skipping array.");

    _mesh.reset(new MeshLib::Mesh(_base_name, nodes, elements, props));
    return true;
}

bool SwmmInterface::matchSubcatchmentsWithPolygons(std::vector<GeoLib::Polyline*> const& lines, std::vector<std::string> const& names)
{
    std::size_t const n_lines (lines.size());
    std::size_t const n_subcatchments (_subcatchments.size());

    if (n_lines != n_subcatchments)
    {
        ERR ("Number of subcatchments does not match number of outlines.");
        return false;
    }
    for (std::size_t i=0; i<n_lines; ++i)
    {
        bool found = false;
        for (std::size_t j=0; j<n_subcatchments; ++j)
        {
            if (names[i] == _subcatchments[j].name)
            {
                _subcatchments[j].outline = lines[i];
                found = true;
                break;
            }
        }
        if (found == false)
        {
            ERR ("No match in subcatcments for outline \"%s\".", names[i].c_str());
            return false;
        }
    }
    return true;
}

std::vector<std::string> SwmmInterface::getSubcatchmentNameMap() const
{
    std::vector<std::string> names;
    names.reserve(_subcatchments.size());
    for (auto sc : _subcatchments)
        names.push_back(sc.name);
    return names;
}

std::vector<std::string> SwmmInterface::getNames(SwmmObject obj_type) const
{
    switch (obj_type)
    {
    case SwmmObject::NODE:
        return _id_nodename_map;
    case SwmmObject::LINK:
        return _id_linkname_map;
    case SwmmObject::SUBCATCHMENT:
        return getSubcatchmentNameMap();
    case SwmmObject::SYSTEM:
        std::vector<std::string> system_name { "System" };
        return system_name;
    }
    ERR ("Object type has no name map");
    std::vector<std::string> empty_vec;
    return empty_vec;
}

std::string SwmmInterface::getName(SwmmObject obj_type, std::size_t idx) const
{
    switch (obj_type)
    {
    case SwmmObject::NODE:
        if (idx < _id_nodename_map.size())
            return _id_nodename_map[idx];
    case SwmmObject::LINK:
        if (idx < _id_linkname_map.size())
            return _id_linkname_map[idx];
    case SwmmObject::SUBCATCHMENT:
        if (idx < _subcatchments.size())
            return _subcatchments[idx].name;
    case SwmmObject::SYSTEM:
        if (idx == 0)
            return std::string("System");
    }
    ERR ("Index out of bounds.");
    return std::string("");
}

std::size_t SwmmInterface::getNumberOfObjects(SwmmObject obj_type) const
{
    std::string const outfile (_base_name + ".out");
    if (OpenSwmmOutFile(const_cast<char*>(outfile.c_str())) != 0)
        return 0;

    switch (obj_type)
    {
        case SwmmObject::SUBCATCHMENT:
            return SWMM_Nsubcatch;
        case SwmmObject::NODE:
            return SWMM_Nnodes;
        case SwmmObject::LINK:
            return SWMM_Nlinks;
        case SwmmObject::SYSTEM:
            return 1;
        default:
            ERR ("Object type not recognised.");
    }
    CloseSwmmOutFile();
    return 0;
}

std::size_t SwmmInterface::getNumberOfParameters(SwmmObject obj_type) const
{
    std::string const outfile (_base_name + ".out");
    std::size_t n_time_steps (std::numeric_limits<std::size_t>::max());
    if (OpenSwmmOutFile(const_cast<char*>(outfile.c_str())) != 0)
        return 0;

    switch (obj_type)
    {
        case SwmmObject::SUBCATCHMENT:
            return (n_obj_params[0] - 1 + SWMM_Npolluts);
        case SwmmObject::NODE:
            return (n_obj_params[1] - 1 + SWMM_Npolluts);
        case SwmmObject::LINK:
            return (n_obj_params[2] - 1 + SWMM_Npolluts);
        case SwmmObject::SYSTEM:
            return n_obj_params[3];
        default:
            ERR ("Object type not recognised.");
    }
    CloseSwmmOutFile();
    return 0;
}

std::size_t SwmmInterface::getNumberOfTimeSteps() const
{
    std::string const outfile (_base_name + ".out");
    if (OpenSwmmOutFile(const_cast<char*>(outfile.c_str())) != 0)
        return std::numeric_limits<std::size_t>::max();
    std::size_t const n_time_steps (static_cast<std::size_t>(SWMM_Nperiods));
    CloseSwmmOutFile();
    return n_time_steps;
}

bool SwmmInterface::addResultsToMesh(MeshLib::Mesh &mesh, SwmmObject const swmm_type,
    std::string const& vec_name, std::vector<double> const& data)
{
    if (!(swmm_type == SwmmObject::NODE) || (swmm_type == SwmmObject::LINK))
    {
        ERR ("Information of this object type cannot be added to mesh.");
        return false;
    }

    if (data.empty())
    {
        ERR ("Data array is empty and cannot be added to mesh.");
        return false;
    }

    if (swmm_type == SwmmObject::NODE && data.size() != mesh.getNumberOfNodes())
    {
        ERR ("Number of mesh nodes (%d) does not match length of array (%d).", mesh.getNumberOfNodes(), data.size());
        return false;
    }

    if (swmm_type == SwmmObject::LINK && data.size() != mesh.getNumberOfElements())
    {
        ERR ("Number of mesh elements (%d) does not match length of array (%d).", mesh.getNumberOfElements(), data.size());
        return false;
    }

    MeshLib::Properties& p = mesh.getProperties();
    MeshLib::MeshItemType item_type = (swmm_type == SwmmObject::NODE) ?
        MeshLib::MeshItemType::Node : MeshLib::MeshItemType::Cell;
    boost::optional<MeshLib::PropertyVector<double>&> prop =
        p.createNewPropertyVector<double>(vec_name, item_type, 1);
    if (!prop)
    {
        ERR ("Error creating array \"%s\".", vec_name.c_str());
        return false;
    }
    prop->reserve(data.size());
    std::copy(data.cbegin(), data.cend(), std::back_inserter(*prop));
    return true;
}

std::vector<double> SwmmInterface::getArrayAtTimeStep(SwmmObject obj_type, std::size_t time_step, std::size_t var_idx) const
{
    std::vector<double> data;
    std::string const outfile (_base_name + ".out");
    if (OpenSwmmOutFile(const_cast<char*>(outfile.c_str())) != 0)
        return data;

    if (time_step >= SWMM_Nperiods)
    {
        ERR ("Time step %d not available, file contains only %d periods.", time_step, SWMM_Nperiods);
        return data;
    }

    bool is_var_idx_okay = true;
    int obj_type_id;
    std::size_t n_objects;
    switch (obj_type)
    {
        case SwmmObject::SUBCATCHMENT:
            obj_type_id = 0;
            n_objects = SWMM_Nsubcatch;
            if (var_idx > (n_obj_params[obj_type_id] - 1 + SWMM_Npolluts))
                is_var_idx_okay = false;
            break;
        case SwmmObject::NODE:
            obj_type_id = 1;
            n_objects = SWMM_Nnodes;
            if (var_idx > (n_obj_params[obj_type_id] + SWMM_Npolluts))
                is_var_idx_okay = false;
            break;
        case SwmmObject::LINK:
            obj_type_id = 2;
            n_objects = SWMM_Nlinks;
            if (var_idx > (n_obj_params[obj_type_id] + SWMM_Npolluts))
                is_var_idx_okay = false;
            break;
        case SwmmObject::SYSTEM:
            obj_type_id = 3;
            n_objects = 1;
            if (var_idx > n_obj_params[obj_type_id])
                is_var_idx_okay = false;
            break;
      default:
         ERR ("Object type not recognised.");
         CloseSwmmOutFile();
         return data;
    }

    if (!is_var_idx_okay)
    {
        ERR ("Requested variable does not exist.");
        CloseSwmmOutFile();
        return data;
    }

    INFO ("Fetching \"%s\"-data for time step %d...",
        getArrayName(obj_type, var_idx, SWMM_Npolluts).c_str(), time_step);

    for (std::size_t i=0; i<n_objects; ++i)
    {
        float val;
        GetSwmmResult(obj_type_id, i, var_idx, time_step, &val);
        data.push_back(static_cast<double>(val));
    }

    CloseSwmmOutFile();
    return data;
}

std::vector<double> SwmmInterface::getArrayForObject(SwmmObject obj_type, std::size_t obj_idx, std::size_t var_idx) const
{
    std::vector<double> data;
    std::string const outfile (_base_name + ".out");
    if (OpenSwmmOutFile(const_cast<char*>(outfile.c_str())) != 0)
        return data;

    bool is_var_idx_okay = true;
    bool is_obj_idx_okay = true;
    int obj_type_id;
    switch (obj_type)
    {
        case SwmmObject::SUBCATCHMENT:
            obj_type_id = 0;
            if (obj_idx >= SWMM_Nsubcatch)
                is_obj_idx_okay = false;
            if (var_idx > (n_obj_params[obj_type_id] + SWMM_Npolluts))
                is_var_idx_okay = false;
            break;
        case SwmmObject::NODE:
            obj_type_id = 1;
            if (obj_idx >= SWMM_Nnodes)
                is_obj_idx_okay = false;
            if (var_idx > (n_obj_params[obj_type_id] + SWMM_Npolluts))
                is_var_idx_okay = false;
            break;
        case SwmmObject::LINK:
            obj_type_id = 2;
            if (obj_idx >= SWMM_Nlinks)
                is_obj_idx_okay = false;
            if (var_idx > (n_obj_params[obj_type_id] + SWMM_Npolluts))
                is_var_idx_okay = false;
            break;
        case SwmmObject::SYSTEM:
            obj_type_id = 3;
            if (obj_idx >= 1)
                is_obj_idx_okay = false;
            if (var_idx > n_obj_params[obj_type_id])
                is_var_idx_okay = false;
            break;
      default:
         ERR ("Object type not recognised.");
         CloseSwmmOutFile();
         return data;
    }

    if (!is_obj_idx_okay)
    {
        ERR ("Requested object index does not exist.");
        CloseSwmmOutFile();
        return data;
    }

    if (!is_var_idx_okay)
    {
        ERR ("Requested variable does not exist.");
        CloseSwmmOutFile();
        return data;
    }

    INFO ("Fetching \"%s\"-data...", getArrayName(obj_type, var_idx, SWMM_Npolluts).c_str());
    std::size_t const n_time_steps (static_cast<std::size_t>(SWMM_Nperiods));
    for (std::size_t i=0; i<n_time_steps; ++i)
    {
        float val;
        GetSwmmResult(obj_type_id , obj_idx, var_idx, i, &val);
        data.push_back(static_cast<double>(val));
    }

    CloseSwmmOutFile();
    return data;
}

std::string SwmmInterface::getArrayName(SwmmObject obj_type, std::size_t var_idx) const
{
    std::string const outfile (_base_name + ".out");
    if (OpenSwmmOutFile(const_cast<char*>(outfile.c_str())) != 0)
        return std::string("");

    std::string const name = getArrayName(obj_type, var_idx, SWMM_Npolluts);
    CloseSwmmOutFile();
    return name;
}

std::string SwmmInterface::getArrayName(SwmmObject obj_type, std::size_t var_idx, std::size_t n_pollutants) const
{
    std::size_t const n_vars (0);
    if (obj_type == SwmmObject::SUBCATCHMENT)
    {
        if (var_idx < n_obj_params[0])
            return subcatchment_vars[var_idx];
        if (var_idx < n_obj_params[0]+n_pollutants)
            return _pollutant_names[var_idx-n_obj_params[0]];
    }
    if (obj_type == SwmmObject::NODE)
    {
        if (var_idx < n_obj_params[1])
            return node_vars[var_idx];
        if (var_idx < n_obj_params[1]+n_pollutants)
            return _pollutant_names[var_idx-n_obj_params[1]];
    }
    if (obj_type == SwmmObject::LINK)
    {
        if (var_idx < n_obj_params[2])
            return link_vars[var_idx];
        if (var_idx < n_obj_params[2]+n_pollutants)
            return _pollutant_names[var_idx-n_obj_params[2]];
    }
    if (obj_type == SwmmObject::SYSTEM && var_idx < n_obj_params[3])
    {
        return system_vars[var_idx];
    }
    ERR ("SwmmInterface::getArrayName() - Index error, no name found.");
    return std::string("");
}

bool SwmmInterface::addRainGaugeTimeSeriesLocations(std::ifstream &in)
{
    std::string line;
    while (getline(in, line))
    {
        if (isSectionFinished(line))
            break;

        if (isCommentLine(line))
            continue;

        std::vector<std::string> const split_str (BaseLib::splitString(line));
        if (split_str.size() < 6)
        {
            ERR ("Rain gauge parameter format not recognised.");
            return false;
        }

        for (auto& stn : _rain_gauges)
        {
            if (stn.first.getName() == split_str[0] && split_str[4] == "FILE")
                stn.second = split_str[5].substr(1, split_str[5].size()-2);
        }
    }

    for (auto const& stn : _rain_gauges)
        if (stn.second.empty())
            WARN ("No associated time series found for rain gauge \"%s\".", stn.first.getName().c_str());
    return true;
}

bool SwmmInterface::readPollutants(std::ifstream &in)
{
    std::string line;
    while (getline(in, line))
    {
        if (isSectionFinished(line))
            return true;

        if (isCommentLine(line))
            continue;

        std::vector<std::string> split_str (BaseLib::splitString(line));
        if (split_str.size() < 6)
        {
            ERR ("Parameter format for pollutants not recognised.");
            return false;
        }
        _pollutant_names.push_back(split_str[0]);
    }
    return true;
}

bool SwmmInterface::isSectionFinished(std::string const& str)
{
    if (str.empty())
        return true;

    std::size_t const pos_beg = str.find_first_not_of(' ', 0);
    if (pos_beg == str.find_first_of(" \n", pos_beg))
        return true;

    return false;
}

bool SwmmInterface::isCommentLine(std::string const& str)
{
    return (str.compare(str.find_first_not_of(' ', 0),1,";") == 0);
}

bool SwmmInterface::getNodeCoordinateVectors(std::vector<double> &x, std::vector<double> &y, std::vector<double> &z) const
{
    std::vector<MeshLib::Node*> const& nodes (_mesh->getNodes());
    for (MeshLib::Node const*const node : nodes)
    {
        x.push_back((*node)[0]);
        y.push_back((*node)[1]);
        z.push_back((*node)[2]);
    }
    return true;
}

bool SwmmInterface::getLinkPointIds(std::vector<std::size_t> &inlets, std::vector<std::size_t> &outlets) const
{
    std::vector<MeshLib::Element*> const& elements (_mesh->getElements());
    for (MeshLib::Element const*const elem : elements)
    {
        if (elem->getGeomType() != MeshLib::MeshElemType::LINE)
        {
            ERR ("Non line-element found in mesh.");
            return false;
        }
        inlets.push_back(elem->getNodeIndex(0));
        outlets.push_back(elem->getNodeIndex(1));
    }
    return true;
}

bool SwmmInterface::writeCsvForTimestep(std::string const& file_name, SwmmObject obj_type, std::size_t time_step) const
{
    FileIO::CsvInterface csv;
    csv.addIndexVectorForWriting(getNumberOfObjects(obj_type));
    csv.addVectorForWriting("Name", getNames(obj_type));
    std::vector<std::string> const obj_names (getNames(obj_type));
    std::size_t const n_params (getNumberOfParameters(obj_type));
    for (std::size_t i=0; i<n_params; ++i)
    {
        std::vector<double> data = getArrayAtTimeStep(obj_type, time_step, i);
        if (!data.empty())
            csv.addVectorForWriting<double>(getArrayName(obj_type, i), data);
    }
    if (csv.getNArrays() < 2)
    {
        ERR ("No data to write");
        return false;
    }
    csv.writeToFile(file_name);
    return true;
}

bool SwmmInterface::writeCsvForObject(std::string const& file_name, SwmmObject obj_type, std::size_t obj_idx) const
{
    FileIO::CsvInterface csv;
    csv.addIndexVectorForWriting(getNumberOfTimeSteps());
    std::size_t const n_params (getNumberOfParameters(obj_type));
    for (std::size_t i=0; i<n_params; ++i)
    {
        std::vector<double> data = getArrayForObject(obj_type, obj_idx, i);
        if (!data.empty())
            csv.addVectorForWriting<double>(getArrayName(obj_type, i), data);
    }
    if (csv.getNArrays() < 2)
    {
        ERR ("No data to write");
        return false;
    }
    csv.writeToFile(file_name);
    return true;
}

} // namespace FileIO
