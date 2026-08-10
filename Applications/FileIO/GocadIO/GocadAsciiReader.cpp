// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "GocadAsciiReader.h"

#include <algorithm>
#include <fstream>
#include <iterator>
#include <memory>
#include <optional>
#include <range/v3/algorithm/max.hpp>
#include <range/v3/algorithm/transform.hpp>
#include <range/v3/view/zip.hpp>
#include <string_view>

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
            // replace chars that will prevent writing the file
            std::replace(mesh_name.begin(), mesh_name.end(), '/', '-');
            std::replace(mesh_name.begin(), mesh_name.end(), '\\', '-');
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
    while (std::getline(in, line))
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

/// Buffers the values of a single node property column while the nodes of one
/// section are parsed.
struct PropertyBuffer
{
    std::vector<double> values;
    MeshLib::PropertyVector<double>* property;
};

/// Parses the node data for the current mesh
bool parseNodes(std::ifstream& in,
                std::vector<MeshLib::Node*>& nodes,
                std::map<std::size_t, std::size_t>& node_id_map,
                MeshLib::Properties const& mesh_prop)
{
    // The buffers follow the iteration over mesh_prop, which is ordered by
    // property name, while the values of a PVRTX line follow the order of the
    // PROPERTIES declaration. Zipping the two below therefore assumes that
    // declaration to be in alphabetical order.
    std::vector<PropertyBuffer> property_buffers;
    for (auto const& [name, property] : mesh_prop)
    {
        if (name == mat_id_name)
        {
            continue;
        }
        if (auto* const p =
                dynamic_cast<MeshLib::PropertyVector<double>*>(property))
        {
            property_buffers.push_back({{}, p});
        }
    }
    // Appends the buffered values of this section to the values of the already
    // parsed sections with a single bulk write per property.
    auto const append_property_buffers = [&property_buffers]()
    {
        for (auto& buffer : property_buffers)
        {
            MeshLib::appendValues(*buffer.property, buffer.values);
        }
    };

    NodeType t = NodeType::UNSPECIFIED;
    std::streampos pos = in.tellg();
    std::string line;
    while (std::getline(in, line))
    {
        std::vector<std::string> str = BaseLib::splitString(line);
        if (line.substr(0, 3) == "SEG" || line.substr(0, 4) == "TRGL")
        {
            in.seekg(pos);
            append_property_buffers();
            return true;
        }

        if (line.substr(0, 28) == "BEGIN_ATOM_REGION_INDICATORS")
        {
            if (!parseAtomRegionIndicators(in))
            {
                ERR("File ended while parsing Atom Region Indicators...");
                return false;
            }
            append_property_buffers();
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
            for (auto& buffer : property_buffers)
            {
                // property_buffers only holds PropertyVector<double> entries,
                // so reading a double per column is safe here.
                double value;
                if (!(sstr >> value))
                {
                    ERR("Error: Could not read the value of property '{:s}' "
                        "of PVRTX line: {:s}",
                        buffer.property->getPropertyName(), line);
                    return false;
                }
                buffer.values.push_back(value);
            }
        }
        else if (line.substr(0, 4) == "ATOM")
        {
            std::size_t new_id;
            std::size_t ref_id;
            std::string keyword;
            sstr >> keyword >> new_id >> ref_id;
            nodes.push_back(new MeshLib::Node(nodes[ref_id]->data(), new_id));
        }
        node_id_map[nodes.back()->getID()] = nodes.size() - 1;
        pos = in.tellg();
    }
    ERR("{:s}", eof_error);
    return false;
}

/// Parses the node ID tuples of all consecutive lines starting with \c keyword,
/// e.g. the segments of a line ("SEG") or the triangles of a surface ("TRGL").
/// The stream is rewound to the first line that does not start with \c keyword.
/// Returns std::nullopt if the file ended while parsing, or if a line does not
/// carry the expected number of node IDs.
template <std::size_t N>
std::optional<std::vector<std::array<std::size_t, N>>> parseNodeIdTuples(
    std::ifstream& in, std::string_view const keyword)
{
    std::vector<std::array<std::size_t, N>> tuples;
    std::streampos pos = in.tellg();
    std::string line;
    while (std::getline(in, line))
    {
        if (line.empty() || isCommentLine(line))
        {
            continue;
        }
        if (!line.starts_with(keyword))
        {
            in.seekg(pos);
            return tuples;
        }
        std::stringstream sstr(line);
        std::string parsed_keyword;
        sstr >> parsed_keyword;
        std::array<std::size_t, N> data{};
        for (auto& node_id : data)
        {
            if (!(sstr >> node_id))
            {
                ERR("Error: Could not read {:d} node IDs of {:s} line: {:s}", N,
                    keyword, line);
                return std::nullopt;
            }
        }
        tuples.push_back(data);
        pos = in.tellg();
    }
    ERR("{:s}", eof_error);
    return std::nullopt;
}

/// Appends \c count material IDs to \c mat_ids for the elements of the Gocad
/// surface or line section that has just been parsed.
///
/// All Gocad sections of a data set share one MaterialIDs vector. Each section
/// gets a fresh material ID, one larger than the current maximum in \c mat_ids
/// (or 0 if \c mat_ids is still empty). The IDs already assigned to earlier
/// sections are not modified.
void extendMaterialIDs(MeshLib::PropertyVector<int>& mat_ids,
                       std::size_t const count)
{
    int const current_mat_id = mat_ids.empty() ? 0 : ranges::max(mat_ids) + 1;
    mat_ids.resize(mat_ids.size() + count, current_mat_id);
}

/// Creates elements from parsed node ID tuples. All elements are constructed
/// first and appended to elems only if every element could be created, i.e.
/// elems is left unchanged on error. The number of nodes per element is taken
/// from the element type.
template <typename ElementType>
bool createElements(
    std::vector<std::array<std::size_t, ElementType::n_all_nodes>> const&
        element_data,
    std::vector<MeshLib::Node*> const& nodes,
    std::vector<MeshLib::Element*>& elems,
    std::map<std::size_t, std::size_t> const& node_id_map)
{
    std::size_t id = elems.size();
    std::vector<std::unique_ptr<ElementType>> new_elems;
    new_elems.reserve(element_data.size());
    for (auto const& data : element_data)
    {
        std::array<MeshLib::Node*, ElementType::n_all_nodes> elem_nodes{};
        for (auto&& [node_id, elem_node] : ranges::views::zip(data, elem_nodes))
        {
            auto const it = node_id_map.find(node_id);
            if (it == node_id_map.end() || it->second >= nodes.size())
            {
                ERR("Error: Node ID ({:d}) out of range [0, {:d}).", node_id,
                    nodes.size());
                return false;
            }
            elem_node = nodes[it->second];
        }
        new_elems.push_back(std::make_unique<ElementType>(elem_nodes, id++));
    }

    elems.reserve(elems.size() + new_elems.size());
    ranges::transform(new_elems, std::back_inserter(elems),
                      [](auto& new_elem) { return new_elem.release(); });

    return true;
}

/// Parses the node ID tuples of all lines starting with \c keyword and creates
/// the corresponding elements including their material IDs. Neither elems nor
/// the material IDs are changed if the elements cannot be created.
template <typename ElementType>
bool parseAndCreateElements(
    std::ifstream& in,
    std::string_view const keyword,
    std::vector<MeshLib::Node*> const& nodes,
    std::vector<MeshLib::Element*>& elems,
    std::map<std::size_t, std::size_t> const& node_id_map,
    MeshLib::Properties& mesh_prop)
{
    // Looked up before the elements are created: a failure here must not
    // leave elems extended without the matching material IDs.
    auto* const mat_ids = mesh_prop.getPropertyVector<int>(mat_id_name);
    if (mat_ids == nullptr)
    {
        ERR("GocadAsciiReader: Property vector '{:s}' not found.", mat_id_name);
        return false;
    }

    auto const element_data =
        parseNodeIdTuples<ElementType::n_all_nodes>(in, keyword);
    if (!element_data)
    {
        return false;
    }
    if (!createElements<ElementType>(*element_data, nodes, elems, node_id_map))
    {
        return false;
    }
    extendMaterialIDs(*mat_ids, element_data->size());
    return true;
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
    if (!parseAndCreateElements<MeshLib::Line>(in, "SEG", nodes, elems,
                                               node_id_map, mesh_prop))
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
    if (!parseAndCreateElements<MeshLib::Tri>(in, "TRGL", nodes, elems,
                                              node_id_map, mesh_prop))
    {
        return false;
    }

    std::string line;
    while (std::getline(in, line))
    {
        std::vector<std::string> str = BaseLib::splitString(line);
        if (str[0] == "TFACE" || str[0] == "3DFace")
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
                          MeshLib::Properties& mesh_prop, T parser,
                          bool const flip_elevation)
{
    std::vector<MeshLib::Node*> nodes;
    std::vector<MeshLib::Element*> elems;
    std::map<std::size_t, std::size_t> node_id_map;
    INFO("Parsing {:s} {:s}.", dataType2ShortString(type), mesh_name);
    bool return_val;
    return_val = parser(in, nodes, elems, node_id_map, mesh_prop);

    if (return_val)
    {
        if (flip_elevation)
        {
            std::for_each(nodes.begin(), nodes.end(),
                          [](MeshLib::Node* n) { (*n)[2] *= -1; });
        }
        return new MeshLib::Mesh(mesh_name, nodes, elems,
                                 true /* compute_element_neighbors */,
                                 mesh_prop);
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
    bool flip_elevation = false;
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
            CoordinateSystem coordinate_system;
            if (!coordinate_system.parse(in))
            {
                ERR("Error parsing coordinate system.");
                return nullptr;
            }
            flip_elevation = (coordinate_system.z_positive ==
                              CoordinateSystem::ZPOSITIVE::Depth);
        }
        else if (str[0] == "GEOLOGICAL_FEATURE" ||
                 str[0] == "GEOLOGICAL_TYPE" ||
                 str[0] == "STRATIGRAPHIC_POSITION" || str[0] == "REGION")
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
            return createMesh(in, type, mesh_name, mesh_prop, parseLine,
                              flip_elevation);
        }
        else if (type == DataType::TSURF &&
                 (str[0] == "TFACE" || str[0] == "3DFace"))
        {
            return createMesh(in, type, mesh_name, mesh_prop, parseSurface,
                              flip_elevation);
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
