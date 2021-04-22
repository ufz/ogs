/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <set>

#include "ProcessLib/ProcessVariable.h"
#include "SecondaryVariable.h"

namespace ProcessLib
{
struct IntegrationPointWriter;
//! Holds information about which variables to write to output files.
struct OutputDataSpecification final
{
    //! All variables that shall be output.
    std::set<std::string> output_variables;

    //! Tells if also to output extrapolation residuals.
    bool const output_residuals;
};

///
/// Prepare the output data, i.e. add the solution to vtu data structure.
void addProcessDataToMesh(
    const double t, std::vector<GlobalVector*> const& x, int const process_id,
    MeshLib::Mesh& mesh,
    std::vector<NumLib::LocalToGlobalIndexMap const*> const& bulk_dof_tables,
    std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_tables,
    std::vector<std::reference_wrapper<ProcessVariable>> const&
        process_variables,
    SecondaryVariableCollection const& secondary_variables,
    bool const output_secondary_variable,
    std::vector<std::unique_ptr<IntegrationPointWriter>> const*
        integration_point_writer,
    OutputDataSpecification const& output_data_specification);

//! Writes output to the given \c file_name using the specified file format.
///
/// See Output::_output_file_data_mode documentation for the data_mode
/// parameter.
enum class OutputType : uint8_t
{
    vtk,
    xdmf
};
void makeOutput(std::string const& file_name, MeshLib::Mesh const& mesh,
                bool const compress_output, int const data_mode,
                OutputType const file_type, int const timestep, double const t);
}  // namespace ProcessLib
