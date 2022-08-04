/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateOutput.h"

#include <memory>
#include <range/v3/algorithm/find.hpp>
#include <tuple>

#include "BaseLib/Algorithm.h"
#include "BaseLib/ConfigTree.h"
#include "BaseLib/FileTools.h"
#include "BaseLib/Logging.h"
#include "BaseLib/cpp23.h"
#include "MeshLib/Mesh.h"
#include "Output.h"

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
}  // namespace

namespace ProcessLib
{
enum class OutputType : uint8_t
{
    vtk,
    xdmf
};

std::unique_ptr<OutputFile> createOutputFile(
    std::string const& output_directory, OutputType const output_type,
    std::string prefix, std::string suffix, std::string const& data_mode,
    bool const compress_output, unsigned int const number_of_files)
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
                compress_output, number_of_files);
        default:
            OGS_FATAL(
                "No supported file type provided. Read '{}' from "
                "<output><type> in prj file. Supported: VTK, XDMF.",
                BaseLib::to_underlying(output_type));
    }
}

std::vector<Output> createOutput(
    const BaseLib::ConfigTree& config,
    std::string const& output_directory,
    std::vector<std::unique_ptr<MeshLib::Mesh>> const& meshes)
{
    DBUG("Parse output configuration:");
    OutputType const output_type = [](auto output_type)
    {
        try
        {
            const std::map<std::string, OutputType> outputType_to_enum = {
                {"VTK", OutputType::vtk}, {"XDMF", OutputType::xdmf}};
            auto type = outputType_to_enum.at(output_type);

            return type;
        }
        catch (std::out_of_range&)
        {
            OGS_FATAL(
                "No supported file type provided. Read `{:s}' from <output><type> \
                in prj File. Supported: VTK, XDMF.",
                output_type);
        }
        //! \ogs_file_param{prj__time_loop__output__type}
    }(config.getConfigParameter<std::string>("type"));

    auto const prefix =
        //! \ogs_file_param{prj__time_loop__output__prefix}
        config.getConfigParameter<std::string>("prefix", "{:meshname}");

    auto const suffix =
        //! \ogs_file_param{prj__time_loop__output__suffix}
        config.getConfigParameter<std::string>("suffix",
                                               "_ts_{:timestep}_t_{:time}");

    auto const compress_output =
        //! \ogs_file_param{prj__time_loop__output__compress_output}
        config.getConfigParameter("compress_output", true);

    auto const hdf =
        //! \ogs_file_param{prj__time_loop__output__hdf}
        config.getConfigSubtreeOptional("hdf");

    auto number_of_files = [&hdf]() -> unsigned int
    {
        if (hdf)
        {
            //! \ogs_file_param{prj__time_loop__output__hdf__number_of_files}
            return hdf->getConfigParameter<unsigned int>("number_of_files");
        }
        return 1;
    }();

    auto const data_mode =
        //! \ogs_file_param{prj__time_loop__output__data_mode}
        config.getConfigParameter<std::string>("data_mode", "Appended");

    // Construction of output times
    std::vector<PairRepeatEachSteps> repeats_each_steps;

    //! \ogs_file_param{prj__time_loop__output__timesteps}
    if (auto const timesteps = config.getConfigSubtreeOptional("timesteps"))
    {
        //! \ogs_file_param{prj__time_loop__output__timesteps__pair}
        for (auto pair : timesteps->getConfigSubtreeList("pair"))
        {
            //! \ogs_file_param{prj__time_loop__output__timesteps__pair__repeat}
            auto repeat = pair.getConfigParameter<unsigned>("repeat");
            //! \ogs_file_param{prj__time_loop__output__timesteps__pair__each_steps}
            auto each_steps = pair.getConfigParameter<unsigned>("each_steps");

            assert(repeat != 0 && each_steps != 0);
            repeats_each_steps.emplace_back(repeat, each_steps);
        }

        if (repeats_each_steps.empty())
        {
            OGS_FATAL(
                "You have not given any pair (<repeat/>, <each_steps/>) that "
                "defines at which timesteps output shall be written. "
                "Aborting.");
        }
    }
    else
    {
        repeats_each_steps.emplace_back(1, 1);
    }

    //! \ogs_file_param{prj__time_loop__output__variables}
    auto const out_vars = config.getConfigSubtree("variables");

    std::set<std::string> output_variables;
    for (auto out_var :
         //! \ogs_file_param{prj__time_loop__output__variables__variable}
         out_vars.getConfigParameterList<std::string>("variable"))
    {
        if (output_variables.find(out_var) != output_variables.cend())
        {
            OGS_FATAL("output variable `{:s}' specified more than once.",
                      out_var);
        }

        DBUG("adding output variable `{:s}'", out_var);
        output_variables.insert(out_var);
    }

    //! \ogs_file_param{prj__time_loop__output__output_extrapolation_residuals}
    bool const output_residuals = config.getConfigParameter<bool>(
        "output_extrapolation_residuals", false);

    std::vector<std::string> mesh_names_for_output;
    //! \ogs_file_param{prj__time_loop__output__meshes}
    if (auto const meshes_config = config.getConfigSubtreeOptional("meshes"))
    {
        if (prefix.find("{:meshname}") == std::string::npos)
        {
            OGS_FATAL(
                "There are multiple meshes defined in the output section of "
                "the project file, but the prefix doesn't contain "
                "'{{:meshname}}'. Thus the names for the files, the simulation "
                "results should be written to, would not be distinguishable "
                "for different meshes.");
        }
        //! \ogs_file_param{prj__time_loop__output__meshes__mesh}
        for (auto mesh_config : meshes_config->getConfigParameterList("mesh"))
        {
            // TODO CL check here that meshes name can be found in meshes list?
            mesh_names_for_output.push_back(
                mesh_config.getValue<std::string>());
            INFO("Configure mesh '{:s}' for output.",
                 mesh_names_for_output.back());
        }
    }

    std::vector<double> fixed_output_times =
        //! \ogs_file_param{prj__time_loop__output__fixed_output_times}
        config.getConfigParameter<std::vector<double>>("fixed_output_times",
                                                       {});
    OutputDataSpecification output_data_specification{
        std::move(output_variables), std::move(fixed_output_times),
        std::move(repeats_each_steps), output_residuals};

    // Remove possible duplicated elements and sort.
    BaseLib::makeVectorUnique(fixed_output_times);

    bool const output_iteration_results =
        //! \ogs_file_param{prj__time_loop__output__output_iteration_results}
        config.getConfigParameter<bool>("output_iteration_results", false);

    auto output_file = createOutputFile(
        output_directory, output_type, std::move(prefix), std::move(suffix),
        data_mode, compress_output, number_of_files);

    std::vector<Output> result;
    result.emplace_back(std::move(output_file), output_iteration_results,
                        std::move(output_data_specification),
                        std::move(mesh_names_for_output), meshes);
    return result;
}

}  // namespace ProcessLib
