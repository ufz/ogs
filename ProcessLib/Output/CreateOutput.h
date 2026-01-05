// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "ProcessLib/Output/Output.h"
#include "ProcessLib/Output/OutputConfig.h"

namespace BaseLib
{
class ConfigTree;
}

namespace MeshLib
{
class Mesh;
}

namespace ProcessLib
{
Output createOutput(OutputConfig&& oc, std::string const& output_directory,
                    std::vector<std::unique_ptr<MeshLib::Mesh>> const& meshes);

std::vector<Output> createOutput(
    const BaseLib::ConfigTree& config,
    const std::string& output_directory,
    std::vector<std::unique_ptr<MeshLib::Mesh>>& meshes);

std::vector<Output> createOutputs(
    const BaseLib::ConfigTree& output_configs,
    std::string const& output_directory,
    std::vector<std::unique_ptr<MeshLib::Mesh>>& meshes);
}  // namespace ProcessLib
