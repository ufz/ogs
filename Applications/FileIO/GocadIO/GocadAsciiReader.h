/**
 *
 * @copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#pragma once

#include <iosfwd>
#include <memory>
#include <string>

#include "MeshLib/Properties.h"

namespace MeshLib
{
    class Mesh;
    class Node;
    class Element;
}

namespace FileIO
{
namespace Gocad
{

enum class GOCAD_DATA_TYPE
{
    UNDEFINED,
    VSET,
    PLINE,
    TSURF,
    MODEL3D
};

class GocadAsciiReader final
{
public:
    /**
     * Constructor takes as argument the Gocad .sg text file.
     */
    explicit GocadAsciiReader();

    GocadAsciiReader(GocadAsciiReader&& src) = delete;
    GocadAsciiReader(GocadAsciiReader const& src) = delete;
    GocadAsciiReader& operator=(GocadAsciiReader&& rhs) = delete;
    GocadAsciiReader& operator=(GocadAsciiReader const& rhs) = delete;

    /// Reads the specified file and writes data into internal mesh vector
    bool readFile(std::string const& file_name, std::vector<std::unique_ptr<MeshLib::Mesh>>& meshes);

private:
    /// Reads one mesh contained in the file (there may be more than one!)
    MeshLib::Mesh* readData(std::ifstream& in, GOCAD_DATA_TYPE const& type, std::string& mesh_name);

    /// Checks if the current line is a comment
    bool isCommentLine(std::string const& str) const;

    /// Checks if a TSurf identifier is found at the current stream position.
    GOCAD_DATA_TYPE datasetFound(std::ifstream& in) const;

    /// Parses the HEADER section (everything except the name is ignored right now)
    bool parseHeader(std::ifstream& in, std::string& mesh_name);

    /// Reads PROPERTY_CLASS_HEADER sections of the file.
    /// All this information is currently ignored.
    bool parsePropertyClass(std::ifstream& in) const;

    /// Parses information of node properties.
    /// Only property names and array sizes are currently used.
    bool parseProperties(std::ifstream& in,
                         std::vector<std::string> const& names,
                         MeshLib::Properties& mesh_prop);

    /// Parses line information (nodes, segments, properties)
    bool parseLine(std::ifstream& in, std::vector<MeshLib::Node*>& nodes,
                   std::vector<MeshLib::Element*>& elems,
                   std::map<std::size_t, std::size_t>& node_id_map,
                   MeshLib::Properties& mesh_prop);

    /// Parses the surface information (nodes, triangles, properties)
    bool parseSurface(std::ifstream& in, std::vector<MeshLib::Node*>& nodes,
                      std::vector<MeshLib::Element*>& elems,
                      std::map<std::size_t, std::size_t>& node_id_map,
                      MeshLib::Properties& mesh_prop);

    /// Parses the node data for the current mesh
    bool parseNodes(std::ifstream& in, std::vector<MeshLib::Node*>& nodes,
                    std::map<std::size_t, std::size_t>& node_id_map,
                    MeshLib::Properties& mesh_prop);

    /// Parses the segments of a line
    bool parseLineSegments(std::ifstream& in, std::vector<MeshLib::Node*>& nodes,
                           std::vector<MeshLib::Element*>& elems,
                           std::map<std::size_t, std::size_t> const& node_id_map,
                           MeshLib::Properties& mesh_prop);

    /// Parses the element data for the current mesh
    bool parseElements(std::ifstream& in, std::vector<MeshLib::Node*>& nodes,
                       std::vector<MeshLib::Element*>& elems,
                       std::map<std::size_t, std::size_t> const& node_id_map,
                       MeshLib::Properties& mesh_prop);

    /// Parses current section until END-tag is reached
    bool skipToEND(std::ifstream& in) const;

    /// Clears the memory if an error occured
    void clearData(std::vector<MeshLib::Node*>& nodes,
                   std::vector<MeshLib::Element*>& elems);

    enum class NODE_TYPE
    {
        UNSPECIFIED,
        VRTX,
        PVRTX
    };
};  // end class GocadTSurfaceReader

}  // end namespace Gocad
}  // end namespace FileIO
