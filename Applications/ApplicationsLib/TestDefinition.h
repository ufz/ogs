/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

namespace BaseLib
{
class ConfigTree;
}

namespace ApplicationsLib
{
class TestDefinition final
{
public:
    /// Constructs test definition from the config and reference path
    /// essentially constructing the command lines to be run on run() function
    /// call.
    TestDefinition(BaseLib::ConfigTree const& config_tree,
                   std::string const& reference_path,
                   std::string const& output_directory);

    bool runTests() const;
    std::vector<std::string> const& getOutputFiles() const;
    std::size_t numberOfTests() const;

private:
    std::vector<std::string> _command_lines;
    std::vector<std::string> _output_files;
};
}  // namespace ApplicationsLib
