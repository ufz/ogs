// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "BaseLib/ExportSymbol.h"

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

    /// Runs all configured test command lines.
    OGS_EXPORT_SYMBOL bool runTests() const;

    /// Runs only test command lines for the given diff tool.
    OGS_EXPORT_SYMBOL bool runTests(std::string_view diff_tool_name) const;

    /// Runs all test command lines except those for the given diff tool.
    OGS_EXPORT_SYMBOL bool runTestsExcluding(
        std::string_view diff_tool_name) const;

    /// Returns true if at least one test command line uses the given diff tool.
    OGS_EXPORT_SYMBOL bool hasTests(std::string_view diff_tool_name) const;

    /// Returns all output files referenced by the configured test definitions.
    std::vector<std::string> const& getOutputFiles() const;

    /// Returns the number of configured test command lines.
    std::size_t numberOfTests() const;

private:
    struct CommandLine
    {
        std::string diff_tool_name;
        std::string command_line;
    };

    static bool runCommandLines(std::vector<CommandLine> const& command_lines,
                                bool fail_on_empty);

    std::vector<CommandLine> _command_lines;
    std::vector<std::string> _output_files;
};
}  // namespace ApplicationsLib
