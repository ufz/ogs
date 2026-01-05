// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
    MeshLib::Mesh* loadMeshFromFile(const std::string& file_name);

    /// Set mesh for writing.
    void setMesh(const MeshLib::Mesh*  mesh);

protected:
    /// Write mesh to stream.
    bool write() override;

private:
    const MeshLib::Mesh* _mesh{nullptr};

};  /* class */

} // end namespace Legacy
} // end namespace IO
} // end namespace MeshLib
