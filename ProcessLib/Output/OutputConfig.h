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

#include <set>
#include <string>
#include <vector>

#include "ProcessLib/Output/OutputDataSpecification.h"

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
}  // namespace ProcessLib
