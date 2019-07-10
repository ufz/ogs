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
class GocadTSurfaceReader final
{
public:
    /**
     * Constructor takes as argument the Gocad .sg text file.
     * @param fname file name
     */
    explicit GocadTSurfaceReader();

    GocadTSurfaceReader(GocadTSurfaceReader&& src) = delete;
    GocadTSurfaceReader(GocadTSurfaceReader const& src) = delete;
    GocadTSurfaceReader& operator=(GocadTSurfaceReader&& rhs) = delete;
    GocadTSurfaceReader& operator=(GocadTSurfaceReader const& rhs) = delete;

    /// Reads the specified file and writes data into internal mesh vector
    bool readFile(std::string const& file_name, std::vector<std::unique_ptr<MeshLib::Mesh>>& meshes);

private:
    /// Reads one mesh contained in the file (there may be more than one!)
    MeshLib::Mesh* readMesh(std::ifstream& in, std::string& mesh_name);

    /// Checks if the current line is a comment
    bool isCommentLine(std::string const& str) const;

    /// Checks if a TSurf identifier is found at the current stream position.
    bool TSurfaceFound(std::ifstream& in) const;

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

    /// Parses the surface information (nodes, triangles, properties)
    bool parseSurface(std::ifstream& in, std::vector<MeshLib::Node*>& nodes,
                      std::vector<MeshLib::Element*>& elems,
                      std::map<std::size_t, std::size_t>& node_id_map,
                      MeshLib::Properties& mesh_prop);

    /// Parses the node data for the current mesh
    bool parseNodes(std::ifstream& in, std::vector<MeshLib::Node*>& nodes,
                    std::map<std::size_t, std::size_t>& node_id_map,
                    MeshLib::Properties& mesh_prop);

    /// Parses the element data for the current mesh
    bool parseElements(std::ifstream& in, std::vector<MeshLib::Node*>& nodes,
                       std::vector<MeshLib::Element*>& elems,
                       std::map<std::size_t, std::size_t> const& node_id_map,
                       MeshLib::Properties& mesh_prop);

    /// Skips over the Model3d sections of the file, should there be any.
    bool skipModel3d(std::ifstream& in) const;

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
