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

#include <stdexcept>
#include <string>

namespace NumLib
{
struct AssemblyException : public std::runtime_error
{
    AssemblyException(std::string const& reason)
        : std::runtime_error{"Error in process' assembly: " + reason} {};
};
}  // namespace NumLib
