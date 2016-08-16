/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_UTILS_PROCESSUTILS_H
#define PROCESSLIB_UTILS_PROCESSUTILS_H

#include <vector>
#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "ProcessLib/Parameter.h"

namespace ProcessLib
{
class ProcessVariable;

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
    BaseLib::ConfigTree const& process_config,
    std::initializer_list<std::string>
        tag_names);

/// Find a parameter of specific type for a name given in the process
/// configuration under the tag.
/// The parameter must have the specified number of components.
/// In the process config a parameter is referenced by a name. For example it
/// will be looking for a parameter named "K" in the list of parameters
/// when the tag is "hydraulic_conductivity":
/// \code
///     <process>
///         ...
///         <hydraulic_conductivity>K</hydraulic_conductivity>
///     </process>
/// \endcode
/// and return a reference to that parameter. Additionally it checks for the
/// type of the found parameter.
template <typename ParameterDataType>
Parameter<ParameterDataType>& findParameter(
    BaseLib::ConfigTree const& process_config, std::string const& tag,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const num_components)
{
    // Find parameter name in process config.
    //! \ogs_file_special
    auto const name = process_config.getConfigParameter<std::string>(tag);

    // Find corresponding parameter by name.
    auto const parameter_it =
        std::find_if(parameters.cbegin(), parameters.cend(),
                     [&name](std::unique_ptr<ParameterBase> const& p) {
                         return p->name == name;
                     });

    if (parameter_it == parameters.end()) {
        OGS_FATAL(
            "Could not find parameter `%s' in the provided parameters list for "
            "config tag <%s>.",
            name.c_str(), tag.c_str());
    }
    DBUG("Found parameter `%s'.", (*parameter_it)->name.c_str());

    // Check the type correctness of the found parameter.
    auto* const parameter =
        dynamic_cast<Parameter<ParameterDataType>*>(parameter_it->get());
    if (!parameter) {
        OGS_FATAL("The read parameter `%s' is of incompatible type.",
                  name.c_str());
    }

    if (parameter->getNumberOfComponents() != num_components) {
        OGS_FATAL(
            "The read parameter `%s' has the wrong number of components (%lu "
            "instead of %u).",
            name.c_str(), parameter->getNumberOfComponents(), num_components);
    }

    return *parameter;
}
}  // ProcessLib

#endif  // PROCESSLIB_UTILS_PROCESSUTILS_H
