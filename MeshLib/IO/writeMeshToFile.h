// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <vtkXMLWriter.h>

#include <filesystem>
#include <set>
#include <string>

namespace MeshLib
{
class Mesh;
}

namespace MeshLib::IO
{
[[nodiscard]] int writeMeshToFile(
    MeshLib::Mesh const& mesh,
    std::filesystem::path const& file_path,
    std::set<std::string> output_variable_names = {},
    bool const use_compression = false,
    int const data_mode = vtkXMLWriter::Appended);
}  // namespace MeshLib::IO
