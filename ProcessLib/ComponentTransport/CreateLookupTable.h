/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

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
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>> const&
        process_variables);

}  // namespace ComponentTransport
}  // namespace ProcessLib
