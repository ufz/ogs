/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>
#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "ProcessLib/Parameter/Parameter.h"

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

/// Find a parameter of specific type for a given name.
///
/// \tparam ParameterDataType the data type of the parameter
/// \param parameter_name name of the requested parameter
/// \param parameters list of parameters in which it will be searched
/// \param num_components the number of components of the parameters or zero if
/// any number is acceptable
///
/// \see The documentation of the other findParameter() function.
template <typename ParameterDataType>
Parameter<ParameterDataType>& findParameter(
    std::string const& parameter_name,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    int const num_components)
{
    // Find corresponding parameter by name.
    auto const parameter_it = std::find_if(
        parameters.cbegin(), parameters.cend(),
        [&parameter_name](std::unique_ptr<ParameterBase> const& p) {
            return p->name == parameter_name;
        });

    if (parameter_it == parameters.end()) {
        OGS_FATAL(
            "Could not find parameter `%s' in the provided parameters list.",
            parameter_name.c_str());
    }

    DBUG("Found parameter `%s'.", (*parameter_it)->name.c_str());

    // Check the type correctness of the found parameter.
    auto* const parameter =
        dynamic_cast<Parameter<ParameterDataType>*>(parameter_it->get());
    if (!parameter) {
        OGS_FATAL("The read parameter `%s' is of incompatible type.",
                  parameter_name.c_str());
    }

    if (num_components != 0 &&
        parameter->getNumberOfComponents() != num_components)
    {
        OGS_FATAL(
            "The read parameter `%s' has the wrong number of components (%lu "
            "instead of %u).",
            parameter_name.c_str(), parameter->getNumberOfComponents(),
            num_components);
    }

    return *parameter;
}

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
    int const num_components)
{
    // Find parameter name in process config.
    //! \ogs_file_special
    auto const name = process_config.getConfigParameter<std::string>(tag);

    return findParameter<ParameterDataType>(name, parameters, num_components);
}
}  // ProcessLib
