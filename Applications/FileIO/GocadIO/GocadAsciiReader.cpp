/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "GocadAsciiReader.h"

#include <iosfwd>

#include "Applications/FileIO/GocadIO/CoordinateSystem.h"
#include "BaseLib/FileTools.h"
#include "BaseLib/Logging.h"
#include "BaseLib/StringTools.h"
#include "MeshLib/Elements/Line.h"
#include "MeshLib/Elements/Tri.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/Properties.h"

namespace FileIO
{
namespace Gocad
{
namespace GocadAsciiReader
{
enum class NodeType
{
    UNSPECIFIED,
    VRTX,
    PVRTX
};

const std::string mat_id_name = "MaterialIDs";
const std::string eof_error = "Error: Unexpected end of file.";

/// A GoCAD file may contain multiple datasets with the same name. To avoid
/// conflicts when writing meshes, a unique id is appended to the mesh name if
/// another dataset with the same name exists.
void checkMeshNames(std::vector<std::unique_ptr<MeshLib::Mesh>> const& meshes)
{
    std::size_t const n_meshes = meshes.size();
    for (std::size_t i = 0; i < n_meshes; ++i)
    {
        std::string const& name = meshes[i]->getName();
        for (std::size_t j = i + 1; j < n_meshes; ++j)
        {
            if (meshes[j]->getName() == name)
            {
                std::string const id_str = std::to_string(meshes[j]->getID());
                meshes[i]->setName(name + "--importID-" + id_str);
                break;
            }
        }
    }
}

/// Checks if the current line is a comment
bool isCommentLine(std::string const& str)
{
    return (str.substr(0, 1) == "#");
}

/// Parses current section until END-tag is reached
bool skipToEND(std::ifstream& in)
{
    std::string line;
    while (std::getline(in, line))
    {
        if (line == "END")
        {
            return true;
        }
    }
    ERR("{:s}", eof_error);
    return false;
}

/// Checks if current line is a designated keyword for a GoCAD data set
bool isKeyword(DataType const t, std::string const& line)
{
    std::size_t str_length = dataType2String(t).length();
    return (line.substr(0, str_length) == dataType2String(t));
}

/// Checks if a GoCAD data set begins at the current stream position.
DataType datasetFound(std::ifstream& in)
{
    std::string line;
    while (std::getline(in, line))
    {
        if (line.empty() || isCommentLine(line))
        {
            continue;
        }

        if (isKeyword(DataType::VSET, line))
        {
            return DataType::VSET;
        }
        if (isKeyword(DataType::PLINE, line))
        {
            return DataType::PLINE;
        }
        if (isKeyword(DataType::TSURF, line))
        {
            return DataType::TSURF;
        }
        if (isKeyword(DataType::MODEL3D, line))
        {
            return DataType::MODEL3D;
        }
        ERR("No known identifier found...");
        return DataType::UNDEFINED;
    }
    return DataType::UNDEFINED;
}
/// Checks if current line is a designated keyword for a GoCAD data set
void checkLineEndings(std::string const& file_name)
{
#ifndef _WIN32
    std::ifstream in(file_name);
    if (in.is_open())
    {
        std::string line;
        std::getline(in, line);
        if (line.back() == '\r')
        {
            OGS_FATAL(
                "Error in input file: {:s}. The line endings are in windows "
                "format. To read this file under UNIX, transform the input "
                "file to unix style line endings (e.g. dos2unix).",
                file_name);
        }
    }
#endif
}

/// Parses the HEADER section (everything except the name is ignored right now)
bool parseHeader(std::ifstream& in, std::string& mesh_name)
{
    std::string line;
    while (std::getline(in, line))
    {
        if (line.substr(0, 5) == "name:")
        {
            mesh_name = line.substr(5, line.length() - 5);
            BaseLib::trim(mesh_name, ' ');
        }
        else if (line.substr(0, 1) == "}")
        {
            return true;
        }
        // ignore all other header parameters
    }
    ERR("{:s}", eof_error);
    return false;
}

/// Reads PROPERTY_CLASS_HEADER sections of the file.
/// All this information is currently ignored.
bool parsePropertyClass(std::ifstream& in)
{
    std::string line;
    while (std::getline(in, line))
    {
        if (line.substr(0, 1) == "}")
        {
            return true;
        }
    }
    ERR("{:s}", eof_error);
    return false;
}

/// Checks if the current line starts with one of the allowed keywords
std::string propertyCheck(std::string const& string)
{
    std::array<std::string, 7> const property_keywords = {
        {"PROPERTY_CLASSES", "PROP_LEGAL_RANGES", "NO_DATA_VALUES",
         "PROPERTY_KINDS", "PROPERTY_SUBCLASSES", "UNITS", "ESIZES"}};

    std::string const str = BaseLib::splitString(string)[0];
    auto res =
        std::find(property_keywords.begin(), property_keywords.end(), str);
    if (res != property_keywords.end())
    {
        return *res;
    }
    return std::string("");
}

/// Parses information of node properties.
/// Only property names and array sizes are currently used.
bool parseProperties(std::ifstream& in,
                     std::vector<std::string> const& names,
                     MeshLib::Properties& mesh_prop)
{
    // Because properties have no end-tag, the position of the last line is
    // stored, so the stream can be set back if none of the allowed property-
    // related keywords is found.
    std::streampos pos = in.tellg();
    std::string line;
    while (getline(in, line))
    {
        std::string const key = propertyCheck(line);
        // This is the intended way to exit this method:
        // No property-related keyword has been found, so the stream is set
        // back one line and the (unrelated) keyword can be read again in the
        // parent method.
        if (key.empty())
        {
            in.seekg(pos);
            return true;
        }

        // Currently all property parameters except array name and size are
        // ignored.
        if (key == "ESIZES")
        {
            std::vector<std::string> prop_size = BaseLib::splitString(line);

            if (names.size() != prop_size.size())
            {
                ERR("Error: Number of PROPERTY-names ({:d}) does not match "
                    "number of ESIZES ({:d})",
                    names.size(), prop_size.size());
                return false;
            }
            std::size_t const n_names(names.size());
            for (std::size_t i = 1; i < n_names; ++i)
            {
                mesh_prop.createNewPropertyVector<double>(
                    names[i],
                    MeshLib::MeshItemType::Node,
                    BaseLib::str2number<std::size_t>(prop_size[i]));
            }
        }
        // Remember current position in case the properties black ends now.
        pos = in.tellg();
    }
    ERR("{:s}", eof_error);
    return false;
}

MeshLib::Node* createNode(std::stringstream& sstr)
{
    std::string keyword;
    std::size_t id;
    std::array<double, 3> data{};
    sstr >> keyword >> id >> data[0] >> data[1] >> data[2];
    return new MeshLib::Node(data, id);
}

/// parse Atom Region Indicator section for current TFACE
/// (currently the information in this section is ignored)
bool parseAtomRegionIndicators(std::ifstream& in)
{
    std::string line;
    while (std::getline(in, line))
    {
        if (line.substr(0, 26) == "END_ATOM_REGION_INDICATORS")
        {
            return true;
        }
    }
    return false;
}

/// Parses the node data for the current mesh
bool parseNodes(std::ifstream& in,
                std::vector<MeshLib::Node*>& nodes,
                std::map<std::size_t, std::size_t>& node_id_map,
                MeshLib::Properties const& mesh_prop)
{
    NodeType t = NodeType::UNSPECIFIED;
    std::streampos pos = in.tellg();
    std::string line;
    while (std::getline(in, line))
    {
        std::vector<std::string> str = BaseLib::splitString(line);
        if (line.substr(0, 3) == "SEG" || line.substr(0, 4) == "TRGL")
        {
            in.seekg(pos);
            return true;
        }

        if (line.substr(0, 28) == "BEGIN_ATOM_REGION_INDICATORS")
        {
            if (!parseAtomRegionIndicators(in))
            {
                ERR("File ended while parsing Atom Region Indicators...");
                return false;
            }
            return true;
        }

        if (line.empty() || isCommentLine(line))
        {
            continue;
        }
        if (!(line.substr(0, 4) == "VRTX" || line.substr(0, 5) == "PVRTX" ||
              line.substr(0, 4) == "ATOM"))
        {
            WARN("GocadAsciiReader::parseNodes() - Unknown keyword found: {:s}",
                 line);
            continue;
        }

        std::stringstream sstr(line);
        if (line.substr(0, 4) == "VRTX" && t != NodeType::PVRTX)
        {
            t = NodeType::VRTX;
            nodes.push_back(createNode(sstr));
        }
        else if (line.substr(0, 5) == "PVRTX" && t != NodeType::VRTX)
        {
            t = NodeType::PVRTX;
            nodes.push_back(createNode(sstr));
            for (auto [name, property] : mesh_prop)
            {
                if (name == mat_id_name)
                {
                    continue;
                }
                if (auto p = dynamic_cast<MeshLib::PropertyVector<double>*>(
                        property))
                {
                    double value;
                    sstr >> value;
                    p->push_back(value);
                }
            }
        }
        else if (line.substr(0, 4) == "ATOM")
        {
            std::size_t new_id;
            std::size_t ref_id;
            std::string keyword;
            sstr >> keyword >> new_id >> ref_id;
            nodes.push_back(
                new MeshLib::Node(nodes[ref_id]->getCoords(), new_id));
        }
        node_id_map[nodes.back()->getID()] = nodes.size() - 1;
        pos = in.tellg();
    }
    ERR("{:s}", eof_error);
    return false;
}

/// Parses the segments of the current line
bool parseLineSegments(std::ifstream& in,
                       std::vector<MeshLib::Node*>& nodes,
                       std::vector<MeshLib::Element*>& elems,
                       std::map<std::size_t, std::size_t> const& node_id_map,
                       MeshLib::Properties& mesh_prop)
{
    MeshLib::PropertyVector<int>& mat_ids =
        *mesh_prop.getPropertyVector<int>(mat_id_name);
    int current_mat_id(0);
    if (!mat_ids.empty())
    {
        current_mat_id = (*std::max_element(mat_ids.begin(), mat_ids.end()))++;
    }
    std::streampos pos = in.tellg();
    std::size_t id(0);
    std::string line;
    while (std::getline(in, line))
    {
        if (line.empty() || isCommentLine(line))
        {
            continue;
        }
        if (line.substr(0, 3) == "SEG")
        {
            std::stringstream sstr(line);
            std::string keyword;
            std::array<std::size_t, 2> data{};
            sstr >> keyword >> data[0] >> data[1];
            std::array<MeshLib::Node*, 2> elem_nodes{};
            for (std::size_t i = 0; i < 2; ++i)
            {
                auto const it = node_id_map.find(data[i]);
                if (it == node_id_map.end() || it->second >= nodes.size())
                {
                    ERR("Error: Node ID ({:d}) out of range (0, {:d}).",
                        data[i], nodes.back()->getID());
                    return false;
                }
                elem_nodes[i] = nodes[it->second];
            }
            elems.push_back(new MeshLib::Line(elem_nodes, id++));
            mat_ids.push_back(current_mat_id);
        }
        else
        {
            in.seekg(pos);
            return true;
        }
        pos = in.tellg();
    }
    ERR("{:s}", eof_error);
    return false;
}

/// Parses line information (nodes, segments, properties)
bool parseLine(std::ifstream& in,
               std::vector<MeshLib::Node*>& nodes,
               std::vector<MeshLib::Element*>& elems,
               std::map<std::size_t, std::size_t>& node_id_map,
               MeshLib::Properties& mesh_prop)
{
    if (!parseNodes(in, nodes, node_id_map, mesh_prop))
    {
        return false;
    }
    if (!parseLineSegments(in, nodes, elems, node_id_map, mesh_prop))
    {
        return false;
    }

    std::string line;
    while (std::getline(in, line))
    {
        std::vector<std::string> str = BaseLib::splitString(line);
        if (str[0] == "ILINE")
        {
            parseLine(in, nodes, elems, node_id_map, mesh_prop);
            return true;
        }
        if (line == "END")
        {
            return true;
        }
        WARN("GocadAsciiReader::parseLine() - Unknown keyword found: {:s}",
             line);
    }
    ERR("{:s}", eof_error);
    return false;
}

/// Parses the element data for the current mesh
bool parseElements(std::ifstream& in,
                   std::vector<MeshLib::Node*>& nodes,
                   std::vector<MeshLib::Element*>& elems,
                   std::map<std::size_t, std::size_t> const& node_id_map,
                   MeshLib::Properties& mesh_prop)
{
    MeshLib::PropertyVector<int>& mat_ids =
        *mesh_prop.getPropertyVector<int>(mat_id_name);
    int current_mat_id(0);
    if (!mat_ids.empty())
    {
        current_mat_id = (*std::max_element(mat_ids.begin(), mat_ids.end()))++;
    }
    std::streampos pos = in.tellg();
    std::size_t id(0);
    std::string line;
    while (std::getline(in, line))
    {
        if (line.empty() || isCommentLine(line))
        {
            continue;
        }
        if (line.substr(0, 4) == "TRGL")
        {
            std::stringstream sstr(line);
            std::string keyword;
            std::array<std::size_t, 3> data{};
            sstr >> keyword >> data[0] >> data[1] >> data[2];
            std::array<MeshLib::Node*, 3> elem_nodes{};
            for (std::size_t i = 0; i < 3; ++i)
            {
                auto const it = node_id_map.find(data[i]);
                if (it == node_id_map.end() || it->second >= nodes.size())
                {
                    ERR("Error: Node ID ({:d}) out of range (0, {:d}).",
                        data[i], nodes.back()->getID());
                    return false;
                }
                elem_nodes[i] = nodes[it->second];
            }
            elems.push_back(new MeshLib::Tri(elem_nodes, id++));
            mat_ids.push_back(current_mat_id);
        }
        else
        {
            in.seekg(pos);
            return true;
        }
        pos = in.tellg();
    }
    ERR("{:s}", eof_error);
    return false;
}

/// Parses the surface information (nodes, triangles, properties)
bool parseSurface(std::ifstream& in,
                  std::vector<MeshLib::Node*>& nodes,
                  std::vector<MeshLib::Element*>& elems,
                  std::map<std::size_t, std::size_t>& node_id_map,
                  MeshLib::Properties& mesh_prop)
{
    if (!parseNodes(in, nodes, node_id_map, mesh_prop))
    {
        return false;
    }
    if (!parseElements(in, nodes, elems, node_id_map, mesh_prop))
    {
        return false;
    }

    std::string line;
    while (std::getline(in, line))
    {
        std::vector<std::string> str = BaseLib::splitString(line);
        if (str[0] == "TFACE")
        {
            parseSurface(in, nodes, elems, node_id_map, mesh_prop);
            return true;
        }
        if (str[0] == "BSTONE")
        {
            // borderstone definition - currently ignored
        }
        else if (str[0] == "BORDER")
        {
            // border tracking direction - currently ignored
        }
        else if (line == "END")
        {
            return true;
        }
        else
        {
            WARN(
                "GocadAsciiReader::parseSurface() - Unknown keyword found: "
                "{:s}",
                line);
        }
    }
    ERR("{:s}", eof_error);
    return false;
}

/// Converts parsed data into mesh
template <typename T>
MeshLib::Mesh* createMesh(std::ifstream& in, DataType type,
                          std::string& mesh_name,
                          MeshLib::Properties& mesh_prop, T parser)
{
    std::vector<MeshLib::Node*> nodes;
    std::vector<MeshLib::Element*> elems;
    std::map<std::size_t, std::size_t> node_id_map;
    INFO("Parsing {:s} {:s}.", dataType2ShortString(type), mesh_name);
    bool return_val;
    return_val = parser(in, nodes, elems, node_id_map, mesh_prop);

    if (return_val)
    {
        return new MeshLib::Mesh(mesh_name, nodes, elems, mesh_prop);
    }
    ERR("Error parsing {:s} {:s}.", dataType2ShortString(type), mesh_name);
    BaseLib::cleanupVectorElements(nodes, elems);
    return nullptr;
}

/// Reads one mesh contained in the file (there may be more than one!)
MeshLib::Mesh* readData(std::ifstream& in,
                        DataType const& type,
                        std::string& mesh_name)
{
    if (!parseHeader(in, mesh_name))
    {
        return nullptr;
    }

    MeshLib::Properties mesh_prop;
    mesh_prop.createNewPropertyVector<int>(mat_id_name,
                                           MeshLib::MeshItemType::Cell, 1);
    std::string line;
    while (std::getline(in, line))
    {
        std::vector<std::string> str = BaseLib::splitString(line);
        if (line.empty() || isCommentLine(line))
        {
            continue;
        }
        if (str[0] == "GOCAD_ORIGINAL_COORDINATE_SYSTEM")
        {
            Gocad::CoordinateSystem coordinate_system;
            if (!coordinate_system.parse(in))
            {
                ERR("Error parsing coordinate system.");
                return nullptr;
            }
        }
        else if (str[0] == "GEOLOGICAL_FEATURE" ||
                 str[0] == "GEOLOGICAL_TYPE" ||
                 str[0] == "STRATIGRAPHIC_POSITION")
        {
            // geological and stratigraphic information - currently ignored
        }
        else if (str[0] == "PROPERTY_CLASS_HEADER")
        {
            if (!parsePropertyClass(in))
            {
                ERR("Error parsing PROPERTY_CLASS_HEADER.");
                return nullptr;
            }
        }
        else if (str[0] == "PROPERTIES")
        {
            if (!parseProperties(in, str, mesh_prop))
            {
                ERR("Error parsing PROPERTIES");
                return nullptr;
            }
        }
        else if (type == DataType::PLINE && str[0] == "ILINE")
        {
            return createMesh(in, type, mesh_name, mesh_prop, parseLine);
        }
        else if (type == DataType::TSURF && str[0] == "TFACE")
        {
            return createMesh(in, type, mesh_name, mesh_prop, parseSurface);
        }
        else
        {
            WARN("GocadAsciiReader::readData() - Unknown keyword found: {:s}",
                 line);
        }
    }
    ERR("{:s}", eof_error);
    return nullptr;
}

bool readFile(std::string const& file_name,
              std::vector<std::unique_ptr<MeshLib::Mesh>>& meshes,
              DataType const export_type)
{
    std::ifstream in(file_name);
    if (!in.is_open())
    {
        ERR("GocadAsciiReader::readFile(): Could not open file {:s}.",
            file_name);
        return false;
    }

    checkLineEndings(file_name);

    DataType type;
    while ((type = datasetFound(in)) != DataType::UNDEFINED)
    {
        if (export_type != DataType::ALL && type != export_type)
        {
            skipToEND(in);
            continue;
        }

        if (type == DataType::VSET || type == DataType::MODEL3D)
        {
            if (!skipToEND(in))
            {
                ERR("Parsing of type {:s} is not implemented. Skipping "
                    "section.",
                    dataType2String(type));
                return false;
            }
            continue;
        }

        std::string mesh_name = BaseLib::dropFileExtension(file_name) +
                                std::to_string(meshes.size() + 1);
        std::unique_ptr<MeshLib::Mesh> mesh(readData(in, type, mesh_name));
        if (mesh == nullptr)
        {
            ERR("File parsing aborted...");
            return false;
        }
        meshes.push_back(std::move(mesh));
    }
    checkMeshNames(meshes);
    return true;
}

}  // namespace GocadAsciiReader
}  // end namespace Gocad
}  // end namespace FileIO
