/**
 *
 * @copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include "GocadTSurfaceReader.h"

#include <logog/include/logog.hpp>

#include "Applications/FileIO/GocadIO/CoordinateSystem.h"
#include "BaseLib/FileTools.h"
#include "BaseLib/StringTools.h"
#include "MeshLib/Elements/Tri.h"
#include "MeshLib/IO/VtkIO/VtuInterface.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"

namespace FileIO
{
namespace Gocad
{
const std::string mat_id_name = "MaterialIDs";
const std::string eof_error = "Error: Unexpected end of file.";

GocadTSurfaceReader::GocadTSurfaceReader(std::string const& filename)
    : _file_name(filename)
{
}

bool GocadTSurfaceReader::readFile()
{
    std::ifstream in(_file_name.c_str());
    if (!in.is_open())
    {
        ERR("GocadTSurfaceReader::allSfc2Mesh(): Could not open file %s.",
            _file_name.c_str());
        return false;
    }

    while (TSurfaceFound(in))
    {
        if (!readMesh(in))
        {
            ERR("File parsing aborted...")
            return false;
        }
    }
    return true;
}

bool GocadTSurfaceReader::readMesh(std::ifstream& in)
{
    std::string const mesh_cnt = "-" + std::to_string(_mesh_vec.size() + 1);
    std::string mesh_name = BaseLib::dropFileExtension(_file_name) + mesh_cnt;
    if (!parseHeader(in, mesh_name))
        return false;

    MeshLib::Properties mesh_prop;
    mesh_prop.createNewPropertyVector<int>(mat_id_name,
                                           MeshLib::MeshItemType::Cell, 1);
    std::string line("");
    while (std::getline(in, line))
    {
        std::vector<std::string> str = BaseLib::splitString(line);
        if (line.empty() || isCommentLine(line))
            continue;
        else if (str[0] == "GOCAD_ORIGINAL_COORDINATE_SYSTEM")
        {
            Gocad::CoordinateSystem coordinate_system;
            if (!coordinate_system.parse(in))
            {
                ERR("Error parsing coordinate system.");
                return false;
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
                return false;
            }
        }
        else if (str[0] == "PROPERTIES")
        {
            if (!parseProperties(in, str, mesh_prop))
            {
                ERR("Error parsing PROPERTIES");
                return false;
            }
        }
        else if (str[0] == "TFACE")
        {
            std::vector<MeshLib::Node*> nodes;
            std::vector<MeshLib::Element*> elems;
            std::map<std::size_t, std::size_t> node_id_map;
            INFO ("Parsing surface %s", mesh_name.c_str());
            if (!parseSurface(in, nodes, elems, node_id_map, mesh_prop))
            {
                ERR("Error parsing Surface %s.", mesh_name.c_str());
                clearData(nodes, elems);
                return false;
            }
            _mesh_vec.push_back(
                new MeshLib::Mesh(mesh_name, nodes, elems, mesh_prop));
            return true;
        }
        else
        {
            WARN("GocadTSurfaceReader::readMesh() - Unknown keyword found: %s",
                 line.c_str());
        }
    }
    ERR("%s", eof_error.c_str());
    return false;
}

MeshLib::Mesh* GocadTSurfaceReader::getData(std::size_t const idx) const
{
    if (_mesh_vec.empty())
    {
        ERR("Error: No mesh data available.");
        return nullptr;
    }
    if (idx < _mesh_vec.size())
        return _mesh_vec[idx];
    ERR("Error: Mesh index (%d) out of range (0, %d).", idx, _mesh_vec.size());
    return nullptr;
}

std::vector<MeshLib::Mesh*> GocadTSurfaceReader::getData() const
{
    if (_mesh_vec.empty())
    {
        ERR("Error: No mesh data available.");
    }
    return _mesh_vec;
}

std::string GocadTSurfaceReader::getMeshName(std::size_t idx) const
{
    if (_mesh_vec.empty())
    {
        ERR("Error: No mesh data available.");
        return std::string();
    }
    if (idx < _mesh_vec.size())
        return _mesh_vec[idx]->getName();
    ERR("Error: Mesh index (%d) out of range (0, %d).", idx, _mesh_vec.size());
    return std::string();
}

void GocadTSurfaceReader::writeData(std::string const& file_name,
                                    std::size_t const idx,
                                    bool write_binary) const
{
    if (_mesh_vec.empty())
    {
        ERR("Error: No mesh data available.");
        return;
    }
    if (idx < _mesh_vec.size())
    {
        int data_mode = (write_binary) ? 2 : 0;
        bool compressed = (write_binary) ? true : false;
        MeshLib::IO::VtuInterface vtu(_mesh_vec[idx], data_mode, compressed);
        vtu.writeToFile(file_name);
        return;
    }
    ERR("Error: Mesh index (%d) out of range (0, %d).", idx, _mesh_vec.size());
    return;
}

void GocadTSurfaceReader::writeData(std::string const& dir,
                                    bool write_binary) const
{
    if (_mesh_vec.empty())
    {
        ERR("Error: No mesh data available.");
        return;
    }
    std::size_t const n_meshes(_mesh_vec.size());
    for (std::size_t i = 0; i < n_meshes; ++i)
    {
        std::string const delim = (dir.back() == '/') ? "" : "/";
        writeData(dir + delim + _mesh_vec[i]->getName(), i, write_binary);
    }
}

bool GocadTSurfaceReader::TSurfaceFound(std::ifstream& in) const
{
    std::string line("");
    while (std::getline(in, line))
    {
        if (line.empty() || isCommentLine(line))
            continue;
        else if (line.substr(0, 11) == "GOCAD TSurf")
            return true;
        // No idea why this is allowed in a *.ts file.
        // It should be a whole different file type.
        else if (line.substr(0, 13) == "GOCAD Model3d")
        {
            if (!skipModel3d(in))
            {
                ERR("Error parsing Model3d");
                return false;
            }
        }
        else
        {
            ERR("No TSurf-identifier found...");
            return false;
        }
    }
    return false;
}

bool GocadTSurfaceReader::isCommentLine(std::string const& str) const
{
    return (str.substr(0, 1) == "#");
}

bool GocadTSurfaceReader::parseHeader(std::ifstream& in, std::string& mesh_name)
{
    std::string line("");
    while (std::getline(in, line))
    {
        if (line.substr(0, 5) == "name:")
            mesh_name = line.substr(5, line.length() - 5);
        else if (line.substr(0, 1) == "}")
            return true;
        // ignore all other header parameters
    }
    ERR("%s", eof_error.c_str());
    return false;
}

bool GocadTSurfaceReader::parsePropertyClass(std::ifstream& in) const
{
    std::string line("");
    while (std::getline(in, line))
    {
        if (line.substr(0, 1) == "}")
            return true;
    }
    ERR("%s", eof_error.c_str());
    return false;
}

/// Checks if the current line starts with one of the allowed keywords
std::string propertyCheck(std::string const& strng)
{
    std::array<std::string, 7> const property_keywords = {
        {"PROPERTY_CLASSES", "PROP_LEGAL_RANGES", "NO_DATA_VALUES",
         "PROPERTY_KINDS", "PROPERTY_SUBCLASSES", "UNITS", "ESIZES"}};

    bool key_found(false);
    std::vector<std::string> str = BaseLib::splitString(strng);
    for (std::string key : property_keywords)
    {
        if (str[0] == key)
            return key;
    }
    return std::string("");
}

bool GocadTSurfaceReader::parseProperties(std::ifstream& in,
                                          std::vector<std::string> const& names,
                                          MeshLib::Properties& mesh_prop)
{
    // Because properties have no end-tag, the position of the last line is
    // stored, so the stream can be set back if none of the allowed property-
    // related keywords is found.
    std::streampos pos = in.tellg();
    std::string line("");
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
                ERR("Error: Number of PROPERTY-names (%d) does not match "
                    "number of ESIZES (%d)", names.size(), prop_size.size());
                return false;
            }
            std::size_t const n_names (names.size());
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
    ERR("%s", eof_error.c_str());
    return false;
}

bool GocadTSurfaceReader::parseSurface(
    std::ifstream& in,
    std::vector<MeshLib::Node*>& nodes,
    std::vector<MeshLib::Element*>& elems,
    std::map<std::size_t, std::size_t>& node_id_map,
    MeshLib::Properties& mesh_prop)
{
    if (!parseNodes(in, nodes, node_id_map, mesh_prop))
        return false;
    if (!parseElements(in, nodes, elems, node_id_map, mesh_prop))
        return false;

    std::string line("");
    while (std::getline(in, line))
    {
        std::vector<std::string> str = BaseLib::splitString(line);
        if (str[0] == "TFACE")
        {
            parseSurface(in, nodes, elems, node_id_map, mesh_prop);
            return true;
        }
        else if (str[0] == "BSTONE")
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
                "GocadTSurfaceReader::parseSurface() - Unknown keyword found: "
                "%s",
                line.c_str());
        }
    }
    ERR("%s", eof_error.c_str());
    return false;
}

MeshLib::Node* createNode(std::stringstream& sstr)
{
    std::string keyword;
    std::size_t id;
    std::array<double, 3> data;
    sstr >> keyword >> id >> data[0] >> data[1] >> data[2];
    return new MeshLib::Node(data, id);
}

bool GocadTSurfaceReader::parseNodes(
    std::ifstream& in,
    std::vector<MeshLib::Node*>& nodes,
    std::map<std::size_t, std::size_t>& node_id_map,
    MeshLib::Properties& mesh_prop)
{
    NODE_TYPE t = NODE_TYPE::UNSPECIFIED;
    double value;
    std::vector<std::string> const array_names =
        mesh_prop.getPropertyVectorNames();
    std::streampos pos = in.tellg();
    std::string line("");
    while (std::getline(in, line))
    {
        std::vector<std::string> str = BaseLib::splitString(line);
        if (line.substr(0, 4) == "TRGL")
        {
            in.seekg(pos);
            return true;
        }

        if (line.empty() || isCommentLine(line))
            continue;
        if (!(line.substr(0, 4) == "VRTX" || line.substr(0, 5) == "PVRTX" ||
              line.substr(0, 4) == "ATOM"))
        {
            WARN(
                "GocadTSurfaceReader::parseNodes() - Unknown keyword found: %s",
                line.c_str());
            continue;
        }

        std::stringstream sstr(line);
        if (line.substr(0, 4) == "VRTX" && t != NODE_TYPE::PVRTX)
        {
            nodes.push_back(createNode(sstr));
        }
        else if (line.substr(0, 5) == "PVRTX" && t != NODE_TYPE::VRTX)
        {
            nodes.push_back(createNode(sstr));
            for (std::string name : array_names)
            {
                if (name == mat_id_name)
                    continue;
                sstr >> value;
                mesh_prop.getPropertyVector<double>(name)->push_back(value);
            }
        }
        else if (line.substr(0, 4) == "ATOM")
        {
            std::size_t new_id, ref_id;
            std::string keyword;
            sstr >> keyword >> new_id >> ref_id;
            nodes.push_back(new MeshLib::Node(nodes[ref_id]->getCoords(), new_id));
        }
        node_id_map[nodes.back()->getID()] = nodes.size() - 1;
        pos = in.tellg();
    }
    ERR("%s", eof_error.c_str());
    return false;
}

bool GocadTSurfaceReader::parseElements(
    std::ifstream& in,
    std::vector<MeshLib::Node*>& nodes,
    std::vector<MeshLib::Element*>& elems,
    std::map<std::size_t, std::size_t> const& node_id_map,
    MeshLib::Properties& mesh_prop)
{
    std::string keyword;
    std::array<std::size_t, 3> data;
    MeshLib::PropertyVector<int>& mat_ids =
        *mesh_prop.getPropertyVector<int>(mat_id_name);
    int current_mat_id(0);
    if (!mat_ids.empty())
        current_mat_id = (*std::max_element(mat_ids.begin(), mat_ids.end()))++;
    std::streampos pos = in.tellg();
    std::size_t id(0);
    std::string line("");
    while (std::getline(in, line))
    {
        if (line.empty() || isCommentLine(line))
            continue;
        if (line.substr(0, 4) == "TRGL")
        {
            std::stringstream sstr(line);
            sstr >> keyword >> data[0] >> data[1] >> data[2];
            std::array<MeshLib::Node*, 3> elem_nodes;
            for (std::size_t i = 0; i < 3; ++i)
            {
                auto const it = node_id_map.find(data[i]);
                if (it == node_id_map.end() || it->second >= nodes.size())
                {
                    ERR("Error: Node ID (%d) out of range (0, %d).", data[i],
                        nodes.back()->getID());
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
    ERR("%s", eof_error.c_str());
    return false;
}

bool GocadTSurfaceReader::skipModel3d(std::ifstream& in) const
{
    std::string line("");
    while (std::getline(in, line))
    {
        if (line == "END")
            return true;
    }
    ERR("%s", eof_error.c_str());
    return false;
}

void GocadTSurfaceReader::clearData(std::vector<MeshLib::Node*>& nodes,
                                    std::vector<MeshLib::Element*>& elems)
{
    for (MeshLib::Element* e : elems)
        delete e;
    for (MeshLib::Node* n : nodes)
        delete n;
}

}  // end namespace Gocad
}  // end namespace FileIO
