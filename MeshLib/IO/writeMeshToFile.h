// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <filesystem>
#include <set>
#include <string>

namespace MeshLib
{
class Mesh;
}

namespace MeshLib::IO
{
int writeMeshToFile(MeshLib::Mesh const& mesh,
                    std::filesystem::path const& file_path,
                    std::set<std::string> variable_output_names = {});
}  // namespace MeshLib::IO