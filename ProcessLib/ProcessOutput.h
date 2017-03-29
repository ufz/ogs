/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ProcessVariable.h"
#include "SecondaryVariable.h"

namespace ProcessLib
{

//! Holds information about which variables to write to output files.
struct ProcessOutput final
{
    //! Constructs a new instance.
    ProcessOutput(BaseLib::ConfigTree const& output_config);

    //! All variables that shall be output.
    std::set<std::string> output_variables;

    //! Tells if also to output extrapolation residuals.
    bool output_residuals = false;
};


//! Writes output to the given \c file_name using the VTU file format.
void doProcessOutput(std::string const& file_name,
                     bool const compress_output,
                     GlobalVector const& x,
                     MeshLib::Mesh& mesh,
                     NumLib::LocalToGlobalIndexMap const& dof_table,
                     std::vector<std::reference_wrapper<ProcessVariable>> const&
                         process_variables,
                     SecondaryVariableCollection secondary_variables,
                     ProcessOutput const& process_output);

} // ProcessLib
