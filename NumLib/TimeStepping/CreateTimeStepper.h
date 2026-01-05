// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
std::unique_ptr<TimeStepAlgorithm> createTimeStepper(
    BaseLib::ConfigTree const& config,
    std::vector<double> const& fixed_times_for_output);
};  // namespace NumLib
