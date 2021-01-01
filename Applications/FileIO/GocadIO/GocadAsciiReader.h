/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "Applications/FileIO/GocadIO/GocadEnums.h"

namespace MeshLib
{
    class Mesh;
}

namespace FileIO
{
namespace Gocad
{
namespace GocadAsciiReader
{

/// Reads the specified file and writes data into internal mesh vector
bool readFile(std::string const& file_name,
              std::vector<std::unique_ptr<MeshLib::Mesh>>& meshes,
              DataType const export_type = DataType::ALL);

}  // namespace GocadAsciiReader
}  // end namespace Gocad
}  // end namespace FileIO
