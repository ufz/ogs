// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <filesystem>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace ProcessLib
{
class ProcessVariable;

namespace ComponentTransport
{
struct LookupTable;

std::unique_ptr<LookupTable> createLookupTable(
    std::optional<std::string> const tabular_file,
    std::filesystem::path const& project_directory,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>> const&
        process_variables);

}  // namespace ComponentTransport
}  // namespace ProcessLib
