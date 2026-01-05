// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <stdexcept>
#include <string>

namespace NumLib
{
struct AssemblyException : public std::runtime_error
{
    explicit AssemblyException(std::string const& reason)
        : std::runtime_error{"Error in process' assembly: " + reason} {};
};
}  // namespace NumLib
