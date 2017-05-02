/**
 * \file
 * \author Karsten Rink
 * \date   2012-05-08
 * \brief  Definition of the MeshIO class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <string>
#include <vector>

#include "BaseLib/IO/Writer.h"
#include "MeshLib/MeshEnums.h"

namespace MeshLib
{
class Mesh;
class Node;
class Element;
template <typename T>
class PropertyVector;
enum class MeshElemType;
namespace IO
{
namespace Legacy
{
/// Interface for handling mesh files from OGS-5 and below. (*.msh files)
class MeshIO : public BaseLib::IO::Writer
{
public:
    /// Constructor.
    MeshIO();

    ~MeshIO() override = default;

    /// Read mesh from file.
    MeshLib::Mesh* loadMeshFromFile(const std::string& fileName);

    /// Set mesh for writing.
    void setMesh(const MeshLib::Mesh*  mesh);

protected:
    /// Write mesh to stream.
    bool write() override;

private:
    void writeElements(std::vector<MeshLib::Element*> const& ele_vec,
                       MeshLib::PropertyVector<int> const* const material_ids,
                       std::ostream& out) const;
    std::size_t readMaterialID(std::istream & in) const;
    MeshLib::Element* readElement(std::istream& line, const std::vector<MeshLib::Node*> &nodes) const;
    std::string ElemType2StringOutput(const MeshLib::MeshElemType t) const;

    const MeshLib::Mesh* _mesh;

};  /* class */

} // end namespace Legacy
} // end namespace IO
} // end namespace MeshLib
