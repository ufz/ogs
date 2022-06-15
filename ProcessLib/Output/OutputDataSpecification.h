/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <set>
#include <string>

namespace ProcessLib
{

//! Holds information about which variables to write to output files.
struct OutputDataSpecification final
{
    OutputDataSpecification(std::set<std::string>&& output_variables,
                            std::vector<double>&& fixed_output_times,
                            bool const output_residuals)
        : output_variables(std::move(output_variables)),
          fixed_output_times(std::move(fixed_output_times)),
          output_residuals(output_residuals)
    {
        if (!std::is_sorted(cbegin(fixed_output_times),
                            cend(fixed_output_times)))
        {
            OGS_FATAL(
                "Vector of fixed output time steps passed to the "
                "OutputDataSpecification constructor must be sorted");
        }
    }

    //! All variables that shall be output.
    std::set<std::string> const output_variables;

    //! Given times that steps have to reach.
    std::vector<double> const fixed_output_times;

    //! Tells if also to output extrapolation residuals.
    bool const output_residuals;
};

}  // namespace ProcessLib
