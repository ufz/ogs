// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>
#include <string>
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
