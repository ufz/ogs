// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
