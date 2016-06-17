/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_PROCESSOUTPUT_H
#define PROCESSLIB_PROCESSOUTPUT_H

#include "ProcessVariable.h"
#include "SecondaryVariable.h"

namespace ProcessLib
{

//! Holds information about which variables to write to output files.
struct ProcessOutput final
{
    //! Constructs a new instance.
    ProcessOutput(BaseLib::ConfigTree const& output_config,
                  std::vector<std::reference_wrapper<ProcessVariable>> const&
                  process_variables,
                  SecondaryVariableCollection const& secondary_variables);

    //! All variables that shall be output.
    std::set<std::string> output_variables;

    //! Tells if also to output extrapolation residuals.
    bool output_residuals = false;

    bool output_iteration_results = false;
};


//! Writes output to the given \c file_name using the VTU file format.
void doProcessOutput(
        std::string const& file_name,
        GlobalVector const& x,
        MeshLib::Mesh& mesh,
        NumLib::LocalToGlobalIndexMap const& dof_table,
        std::vector<std::reference_wrapper<ProcessVariable>> const&
        process_variables,
        SecondaryVariableCollection secondary_variables,
        ProcessOutput const& process_output);

} // ProcessLib


#endif // PROCESSLIB_PROCESSOUTPUT_H
