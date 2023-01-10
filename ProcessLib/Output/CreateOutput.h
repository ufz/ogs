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

#include "ProcessLib/Output/Output.h"
#include "ProcessLib/Output/OutputDataSpecification.h"

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
enum class OutputType : uint8_t
{
    vtk,
    xdmf
};

struct OutputConfig
{
    OutputType output_type;
    std::string prefix;
    std::string suffix;
    bool compress_output;
    unsigned int number_of_files;
    std::string data_mode;
    std::vector<PairRepeatEachSteps> repeats_each_steps;
    std::set<std::string> output_variables;
    bool output_extrapolation_residuals;
    std::vector<std::string> mesh_names_for_output;
    std::vector<double> fixed_output_times;
    bool output_iteration_results;
};

OutputConfig createOutputConfig(const BaseLib::ConfigTree& config);

Output createOutput(OutputConfig&& oc, std::string const& output_directory,
                    std::vector<std::unique_ptr<MeshLib::Mesh>> const& meshes);

std::vector<Output> createOutput(
    const BaseLib::ConfigTree& config,
    const std::string& output_directory,
    std::vector<std::unique_ptr<MeshLib::Mesh>> const& meshes);

std::vector<Output> createOutputs(
    const BaseLib::ConfigTree& config,
    std::string const& output_directory,
    std::vector<std::unique_ptr<MeshLib::Mesh>> const& meshes);
}  // namespace ProcessLib
