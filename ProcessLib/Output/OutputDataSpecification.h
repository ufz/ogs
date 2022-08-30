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

#include <spdlog/fmt/bundled/ostream.h>

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

inline std::ostream& operator<<(std::ostream& os,
                                PairRepeatEachSteps const& pair)
{
    os << "Output " << pair.repeat << " times every " << pair.each_steps
       << " timestep.\n";
    return os;
}

//! Holds information about which variables to write to output files.
struct OutputDataSpecification final
{
    OutputDataSpecification(
        std::set<std::string>&& output_variables_,
        std::vector<double>&& fixed_output_times_,
        std::vector<PairRepeatEachSteps>&& repeats_each_steps_,
        bool const output_residuals_)
        : output_variables(std::move(output_variables_)),
          fixed_output_times(std::move(fixed_output_times_)),
          repeats_each_steps(std::move(repeats_each_steps_)),
          output_residuals(output_residuals_)
    {
        if (!std::is_sorted(cbegin(fixed_output_times),
                            cend(fixed_output_times)))
        {
            OGS_FATAL(
                "Vector of fixed output time steps passed to the "
                "OutputDataSpecification constructor must be sorted");
        }
        // check the repeats_each_steps pairs
        for (auto const& pair : repeats_each_steps)
        {
            if (pair.each_steps == 0)
            {
                OGS_FATAL(
                    "Step in pair of <repeats><steps> is zero but has to be "
                    "greater than zero.");
            }
        }
        if (repeats_each_steps.empty())
        {
            repeats_each_steps.emplace_back(1, std::numeric_limits<int>::max());
        }
    }

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
    bool isOutputStep(int timestep, double const time) const
    {
        auto isFixedOutputStep = [this](double const time) -> bool
        {
            auto const fixed_output_time = std::lower_bound(
                cbegin(fixed_output_times), cend(fixed_output_times), time);
            return ((fixed_output_time != cend(fixed_output_times)) &&
                    (std::abs(*fixed_output_time - time) <
                     std::numeric_limits<double>::epsilon()));
        };

        auto isPairRepeatsEachTimeStepOutput = [this](int timestep) -> bool
        {
            int each_steps = 1;

            for (auto const& pair : repeats_each_steps)
            {
                each_steps = pair.each_steps;

                if (timestep > pair.repeat * each_steps)
                {
                    timestep -= pair.repeat * each_steps;
                }
                else
                {
                    break;
                }
            }

            return timestep % each_steps == 0;
        };

        return isFixedOutputStep(time) ||
               isPairRepeatsEachTimeStepOutput(timestep);
    }
};

inline std::ostream& operator<<(std::ostream& os,
                                OutputDataSpecification const& o)
{
    os << "OuputDataSpecification" << std::endl;
    os << "\toutput_variables: ";
    std::copy(o.output_variables.begin(), o.output_variables.end(),
              std::ostream_iterator<std::string>(os, " "));
    os << "\n";
    os << "\tfixed_output_times: ";
    std::copy(o.fixed_output_times.begin(), o.fixed_output_times.end(),
              std::ostream_iterator<double>(os, " "));
    os << "\n";
    os << "\trepeats_each_steps: ";
    std::copy(o.repeats_each_steps.begin(), o.repeats_each_steps.end(),
              std::ostream_iterator<PairRepeatEachSteps>(os, " "));
    os << "\n";
    os << "\toutput_residual: " << o.output_residuals << "\n";
    return os;
}

}  // namespace ProcessLib
