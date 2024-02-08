/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on November 21, 2023, 3:37 PM
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
class ConvergenceCriterion;

class StaggeredCoupling;

struct LocalCouplingParameters
{
    std::vector<std::string> process_names;
    int max_iterations;
};

std::tuple<std::vector<std::unique_ptr<NumLib::ConvergenceCriterion>>,
           std::vector<LocalCouplingParameters>,
           int>
parseCoupling(BaseLib::ConfigTree const& config);

/// Create a StaggeredCoupling instance from the given configuration.
template <typename ProcessData>
std::unique_ptr<StaggeredCoupling> createStaggeredCoupling(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ProcessData>> const& per_process_data);

}  // namespace NumLib

#include "CreateStaggeredCoupling-impl.h"
