// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "FEFLOWGeoInterface.h"

#include <libxml/parser.h>
#include <libxml/tree.h>

#include <boost/algorithm/string/trim.hpp>
#include <cctype>
#include <fstream>
#include <memory>
#include <sstream>

#include "BaseLib/FileTools.h"
#include "BaseLib/Logging.h"
#include "BaseLib/StringTools.h"
#include "GeoLib/GEOObjects.h"
#include "GeoLib/Point.h"
#include "GeoLib/Polygon.h"

namespace
{

inline std::string getProp(xmlNodePtr node, const xmlChar* name)
{
    if (!node)
    {
        return {};
    }
    xmlChar* v = xmlGetProp(node, name);
    if (!v)
    {
        return {};
    }
    std::string s(reinterpret_cast<const char*>(v));
    xmlFree(v);
    return s;
}

inline std::string getText(xmlNodePtr node)
{
    if (!node)
    {
        return {};
    }
    xmlChar* v = xmlNodeGetContent(node);
    if (!v)
    {
        return {};
    }
    std::string s(reinterpret_cast<const char*>(v));
    xmlFree(v);
    return s;
}

inline xmlNodePtr firstChild(xmlNodePtr parent, const xmlChar* name)
{
    if (!parent)
    {
        return nullptr;
    }
    for (xmlNodePtr n = parent->children; n; n = n->next)
    {
        if (n->type == XML_ELEMENT_NODE && xmlStrcmp(n->name, name) == 0)
        {
            return n;
        }
    }
    return nullptr;
}
}  // namespace

namespace FileIO
{
void FEFLOWGeoInterface::readFEFLOWFile(const std::string& filename,
                                        GeoLib::GEOObjects& geo_objects)
{
    std::ifstream in(filename.c_str());
    if (!in)
    {
        ERR("FEFLOWGeoInterface::readFEFLOWFile(): Could not open file {:s}.",
            filename);
        return;
    }

    unsigned dimension = 0;
    std::vector<GeoLib::Point*> points;
    std::vector<GeoLib::Polyline*> lines;

    bool isXZplane = false;

    std::string line_string;
    std::stringstream line_stream;
    while (!in.eof())
    {
        std::getline(in, line_string);
        if (line_string.find("CLASS") != std::string::npos)
        {
            std::getline(in, line_string);
            line_stream.str(line_string);
            unsigned dummy = 0;
            for (int i = 0; i < 3; i++)
            {
                line_stream >> dummy;
            }
            line_stream >> dimension;
            line_stream.clear();
        }
        else if (line_string == "GRAVITY")
        {
            std::getline(in, line_string);
            line_stream.str(line_string);
            double vec[3] = {};
            line_stream >> vec[0] >> vec[1] >> vec[2];
            if (vec[0] == 0.0 && vec[1] == -1.0 && vec[2] == 0.0)
            {
                isXZplane = true;
            }
            line_stream.clear();
        }
        else if (line_string == "SUPERMESH")
        {
            readSuperMesh(in, dimension, points, lines);
        }
    }
    in.close();

    if (isXZplane && !points.empty())
    {
        for (auto* pt : points)
        {
            (*pt)[2] = (*pt)[1];
            (*pt)[1] = .0;
        }
    }
    std::string project_name(
        BaseLib::extractBaseNameWithoutExtension(filename));
    if (!points.empty())
    {
        geo_objects.addPointVec(std::move(points), project_name,
                                GeoLib::PointVec::NameIdMap{});
    }
    if (!lines.empty())
    {
        geo_objects.addPolylineVec(std::move(lines), project_name,
                                   GeoLib::PolylineVec::NameIdMap{});
    }
}

void FEFLOWGeoInterface::readPoints(xmlNodePtr nodesEle,
                                    const xmlChar* tag,
                                    int dim,
                                    std::vector<GeoLib::Point*>& points)
{
    xmlNodePtr xmlEle = firstChild(nodesEle, tag);
    if (xmlEle == nullptr)
    {
        return;
    }

    std::istringstream ss(getText(xmlEle));
    std::string line_str;
    while (std::getline(ss, line_str))
    {
        boost::algorithm::trim_right(line_str);
        if (line_str.empty())
        {
            continue;
        }
        std::istringstream line_ss(line_str);
        std::size_t pt_id = 0;
        std::array<double, 3> pt_xyz{};
        line_ss >> pt_id;
        for (int i = 0; i < dim; ++i)
        {
            line_ss >> pt_xyz[i];
        }
        points[pt_id - 1] = new GeoLib::Point(pt_xyz, pt_id);
    }
}

void FEFLOWGeoInterface::readSuperMesh(std::ifstream& in,
                                       unsigned dimension,
                                       std::vector<GeoLib::Point*>& points,
                                       std::vector<GeoLib::Polyline*>& lines)
{
    // Read XML block
    std::ostringstream oss;
    std::string line_string;
    while (true)
    {
        std::getline(in, line_string);
        BaseLib::trim(line_string);
        oss << line_string << "\n";
        if (line_string.find("</supermesh>") != std::string::npos)
        {
            break;
        }
    }
    const std::string xml_str = oss.str();

    xmlInitParser();
    xmlDocPtr doc =
        xmlReadMemory(xml_str.c_str(), static_cast<int>(xml_str.size()),
                      "supermesh.xml", nullptr, XML_PARSE_NOBLANKS);
    if (doc == nullptr)
    {
        ERR("FEFLOWGeoInterface::readSuperMesh(): Illegal XML format error");
        xmlCleanupParser();
        return;
    }

    xmlNodePtr docElem = xmlDocGetRootElement(doc);
    if (docElem == nullptr)
    {
        xmlFreeDoc(doc);
        xmlCleanupParser();
        return;
    }

    static constexpr xmlChar nodes_name[] = "nodes";
    xmlNodePtr nodesEle = firstChild(docElem, nodes_name);
    if (nodesEle == nullptr)
    {
        xmlFreeDoc(doc);
        xmlCleanupParser();
        return;
    }

    static constexpr xmlChar count_name[] = "count";
    const std::string cnt = getProp(nodesEle, count_name);
    const long n_points = cnt.empty() ? 0L : std::stol(cnt);

    // Use unique_ptr for RAII-based cleanup
    std::vector<std::unique_ptr<GeoLib::Point>> temp_points;
    temp_points.resize(static_cast<std::size_t>(n_points));

    // Convert to raw pointer vector for readPoints compatibility
    std::vector<GeoLib::Point*> raw_points(static_cast<std::size_t>(n_points),
                                           nullptr);

    static constexpr xmlChar fixed_name[] = "fixed";
    readPoints(nodesEle, fixed_name, static_cast<int>(dimension), raw_points);
    static constexpr xmlChar linear_name[] = "linear";
    readPoints(nodesEle, linear_name, static_cast<int>(dimension), raw_points);
    static constexpr xmlChar parabolic_name[] = "parabolic";
    readPoints(nodesEle, parabolic_name, static_cast<int>(dimension),
               raw_points);

    // Transfer ownership to unique_ptrs
    for (std::size_t i = 0; i < raw_points.size(); ++i)
    {
        temp_points[i].reset(raw_points[i]);
    }

    static constexpr xmlChar polygons_name[] = "polygons";
    xmlNodePtr polygonsEle = firstChild(docElem, polygons_name);
    if (polygonsEle == nullptr)
    {
        // unique_ptrs will automatically clean up
        xmlFreeDoc(doc);
        xmlCleanupParser();
        return;
    }

    std::vector<std::unique_ptr<GeoLib::Polyline>> temp_lines;

    for (xmlNodePtr child = polygonsEle->children; child; child = child->next)
    {
        static constexpr xmlChar polygon_name[] = "polygon";
        if (child->type != XML_ELEMENT_NODE ||
            xmlStrcmp(child->name, polygon_name) != 0)
        {
            continue;
        }

        xmlNodePtr nodes = firstChild(child, nodes_name);
        if (nodes == nullptr)
        {
            continue;
        }

        const std::string cnt_nodes = getProp(nodes, count_name);
        const std::size_t n_pts =
            cnt_nodes.empty() ? 0u
                              : static_cast<std::size_t>(std::stol(cnt_nodes));
        std::string pt_ids = getText(nodes);
        BaseLib::trim(pt_ids);

        auto line = std::make_unique<GeoLib::Polyline>(raw_points);

        std::istringstream ss(pt_ids);
        for (std::size_t i = 0; i < n_pts; ++i)
        {
            int pt_id = 0;
            ss >> pt_id;
            line->addPoint(pt_id - 1);
        }
        // close the polygon
        line->addPoint(line->getPointID(0));

        temp_lines.push_back(std::move(line));
    }

    xmlFreeDoc(doc);
    xmlCleanupParser();

    // Only transfer ownership to output vectors if everything succeeded
    for (auto& pt : temp_points)
    {
        points.push_back(pt.release());
    }
    for (auto& line : temp_lines)
    {
        lines.push_back(line.release());
    }
}

}  // namespace FileIO
