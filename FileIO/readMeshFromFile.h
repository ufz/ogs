/**
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * @file readMeshFromFile.h
 * @date 2012-09-27
 * @author Karsten Rink
 */

#ifndef READMESHFROMFILE_H
#define READMESHFROMFILE_H

#include <string>

namespace MeshLib
{
class Mesh;
}

namespace FileIO
{
MeshLib::Mesh* readMeshFromFile(const std::string &file_name);
}

#endif // READMESHFROMFILE_H
