/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <ostream>
#include <set>
#include <string>
#include <vector>

#include "BaseLib/Error.h"

namespace ProcessLib
{
struct PairRepeatEachSteps
{
    explicit PairRepeatEachSteps(int c, int e) : repeat(c), each_steps(e) {}

    const int repeat;      //!< Apply \c each_steps \c repeat times.
    const int each_steps;  //!< Do output every \c each_steps timestep.
};

std::ostream& operator<<(std::ostream& os, PairRepeatEachSteps const& pair);

//! Holds information about which variables to write to output files.
struct OutputDataSpecification final
{
    OutputDataSpecification(
        std::set<std::string>&& output_variables_,
        std::vector<double>&& fixed_output_times_,
        std::vector<PairRepeatEachSteps>&& repeats_each_steps_,
        bool const output_residuals_);

    //! All variables that shall be output.
    std::set<std::string> output_variables;

    //! Given times that steps have to reach.
    std::vector<double> fixed_output_times;

    //! Describes after which timesteps to write output.
    std::vector<PairRepeatEachSteps> repeats_each_steps;

    //! Tells if also to output extrapolation residuals.
    bool output_residuals;

    //! Determines if there should be output at the given \c timestep or \c
    //! time.
    bool isOutputStep(int timestep, double const time) const;
};

std::ostream& operator<<(std::ostream& os, OutputDataSpecification const& o);

}  // namespace ProcessLib
