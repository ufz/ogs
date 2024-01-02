/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateOutput.h"

#include <memory>
#include <range/v3/algorithm/find.hpp>
#include <tuple>

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Logging.h"
#include "BaseLib/cpp23.h"
#include "MeshLib/Mesh.h"
#include "ProcessLib/Output/CreateOutputConfig.h"

namespace
{
//! Converts a vtkXMLWriter's data mode string to an int. See
/// OutputVTKFormat::data_mode.
int convertVtkDataMode(std::string_view const& data_mode)
{
    using namespace std::string_view_literals;
    constexpr std::array data_mode_lookup_table{"Ascii"sv, "Binary"sv,
                                                "Appended"sv};
    auto res = ranges::find(begin(data_mode_lookup_table),
                            end(data_mode_lookup_table), data_mode);
    if (res == data_mode_lookup_table.end())
    {
        OGS_FATAL(
            "Unsupported vtk output file data mode '{:s}'. Expected Ascii, "
            "Binary, or Appended.",
            data_mode);
    }
    return static_cast<int>(std::distance(begin(data_mode_lookup_table), res));
}

bool areOutputNamesUnique(std::vector<ProcessLib::Output> const& outputs)
{
    std::vector<std::string> output_names;
    for (auto const& output : outputs)
    {
        auto output_mesh_names = output.getFileNamesForOutput();
        output_names.insert(output_names.end(), output_mesh_names.begin(),
                            output_mesh_names.end());
    }
    std::sort(output_names.begin(), output_names.end());
    auto const last = std::unique(output_names.begin(), output_names.end());
    return last == output_names.end();
}
}  // namespace

namespace ProcessLib
{
std::unique_ptr<OutputFormat> createOutputFormat(
    std::string const& output_directory, OutputType const output_type,
    std::string prefix, std::string suffix, std::string const& data_mode,
    bool const compress_output, unsigned int const number_of_files,
    unsigned int const chunk_size_bytes)
{
    switch (output_type)
    {
        case OutputType::vtk:
            return std::make_unique<OutputVTKFormat>(
                output_directory, std::move(prefix), std::move(suffix),
                compress_output, convertVtkDataMode(data_mode));
        case OutputType::xdmf:
            return std::make_unique<OutputXDMFHDF5Format>(
                output_directory, std::move(prefix), std::move(suffix),
                compress_output, number_of_files, chunk_size_bytes);
        default:
            OGS_FATAL(
                "No supported file type provided. Read '{}' from "
                "<output><type> in prj file. Supported: VTK, XDMF.",
                BaseLib::to_underlying(output_type));
    }
}

Output createOutput(OutputConfig&& oc, std::string const& output_directory,
                    std::vector<std::unique_ptr<MeshLib::Mesh>> const& meshes)
{
    auto output_format = createOutputFormat(
        output_directory, oc.output_type, std::move(oc.prefix),
        std::move(oc.suffix), oc.data_mode, oc.compress_output,
        oc.number_of_files, oc.chunk_size_bytes);

    OutputDataSpecification output_data_specification{
        std::move(oc.output_variables), std::move(oc.fixed_output_times),
        std::move(oc.repeats_each_steps), oc.output_extrapolation_residuals};

    return {std::move(output_format), oc.output_iteration_results,
            std::move(output_data_specification),
            std::move(oc.mesh_names_for_output), meshes};
}

std::vector<Output> createOutput(
    const BaseLib::ConfigTree& config,
    std::string const& output_directory,
    std::vector<std::unique_ptr<MeshLib::Mesh>>& meshes)
{
    std::vector<Output> outputs;
    auto oc = createOutputConfig(config, meshes);
    outputs.push_back(createOutput(std::move(oc), output_directory, meshes));
    return outputs;
}

std::vector<Output> createOutputs(
    const BaseLib::ConfigTree& output_configs,
    std::string const& output_directory,
    std::vector<std::unique_ptr<MeshLib::Mesh>>& meshes)
{
    DBUG("Parse outputs configuration:");
    std::vector<Output> outputs;
    for (auto const& output_config :
         //! \ogs_file_param{prj__time_loop__outputs__output}
         output_configs.getConfigSubtreeList("output"))
    {
        auto oc = createOutputConfig(output_config, meshes);
        outputs.push_back(
            createOutput(std::move(oc), output_directory, meshes));
    }
    if (areOutputNamesUnique(outputs))
    {
        return outputs;
    }
    else
    {
        OGS_FATAL(
            "Output configuration paths are not unique. This will lead to "
            "overwritten results or invalid / corrupted data within the "
            "files.");
    }
}
}  // namespace ProcessLib
