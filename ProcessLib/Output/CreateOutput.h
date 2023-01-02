/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

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
class Process;
class Output;
}  // namespace ProcessLib

namespace ProcessLib
{
std::vector<Output> createOutput(
    const BaseLib::ConfigTree& config,
    const std::string& output_directory,
    std::vector<std::unique_ptr<MeshLib::Mesh>> const& meshes);

std::vector<Output> createOutputs(
    const BaseLib::ConfigTree& config,
    std::string const& output_directory,
    std::vector<std::unique_ptr<MeshLib::Mesh>> const& meshes);
}  // namespace ProcessLib
