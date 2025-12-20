// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <string>

namespace MeshLib
{
class Mesh;
}

namespace FileIO
{
namespace GMSH
{
/**
 * checks if there is a GMSH mesh file header
 * \param fname the file name of the mesh (including the path)
 * \return true, if the file seems to be a valid GMSH file, else false
 */
bool isGMSHMeshFile(const std::string& fname);

/**
 * reads a mesh created by GMSH - this implementation is based on the former
 * function GMSH2MSH
 * \param fname the file name of the mesh (including the path)
 * \param is_created_with_gmsh2 An indicator for the mesh created by using Gmsh
 * version 2.
 * \return
 */
MeshLib::Mesh* readGMSHMesh(std::string const& fname,
                            bool const is_created_with_gmsh2 = false);

}  // end namespace GMSH
}  // end namespace FileIO
