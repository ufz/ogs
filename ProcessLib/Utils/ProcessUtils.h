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

#include <string>
#include <vector>

namespace BaseLib
{
class ConfigTree;
}
namespace ProcessLib
{
class ProcessVariable;
}

namespace ProcessLib
{
/// Find process variables in \c variables whose names match the settings under
/// the given \c tag_names in the \c process_config.
///
/// In the process config a process variable is referenced by a name. For
/// example it will be looking for a variable named "H" in the list of process
/// variables when the tag is "hydraulic_head":
/// \code
///     <process>
///         ...
///         <process_variables>
///             <hydraulic_head>H</hydraulic_head>
///             ...
///         </process_variables>
///         ...
///     </process>
/// \endcode
///
/// \return a vector of references to the found variable(s).
std::vector<std::reference_wrapper<ProcessVariable>> findProcessVariables(
    std::vector<ProcessVariable> const& variables,
    BaseLib::ConfigTree const& pv_config,
    std::initializer_list<std::string>
        tags);

std::vector<std::reference_wrapper<ProcessVariable>> findProcessVariables(
    std::vector<ProcessVariable> const& variables,
    BaseLib::ConfigTree const& pv_config,
    std::string const& tag);

}  // namespace ProcessLib
