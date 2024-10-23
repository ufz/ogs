/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  Created on June 26, 2017, 4:43 PM
 */

#pragma once

#include <memory>
#include <vector>

namespace BaseLib
{
class ConfigTree;
}

namespace NumLib
{
class TimeStepAlgorithm;

struct EvolutionaryPIDcontrollerParameters final
{
    double t0;
    double t_end;
    double h0;
    double h_min;
    double h_max;
    double rel_h_min;
    double rel_h_max;
    double tol;
};

/// Parse an EvolutionaryPIDcontroller time stepper from the given
/// configuration
EvolutionaryPIDcontrollerParameters parseEvolutionaryPIDcontroller(
    BaseLib::ConfigTree const& config);

/// Create an EvolutionaryPIDcontroller time stepper from the given
/// configuration
std::unique_ptr<TimeStepAlgorithm> createEvolutionaryPIDcontroller(
    EvolutionaryPIDcontrollerParameters const& config,
    std::vector<double> const& fixed_times_for_output);
}  // end of namespace NumLib
