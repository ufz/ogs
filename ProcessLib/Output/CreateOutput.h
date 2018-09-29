/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <string>

namespace BaseLib
{
class ConfigTree;
}

namespace ProcessLib
{
class Process;
class Output;
}  // namespace ProcessLib

namespace ProcessLib
{
std::unique_ptr<Output> createOutput(const BaseLib::ConfigTree& config,
                                     const std::string& output_directory);

}  // namespace ProcessLib
