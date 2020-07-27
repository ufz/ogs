/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>
#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"

#include "Parameter.h"

namespace ParameterLib
{
/// Find an optional parameter of specific type for a given name.
///
/// \tparam ParameterDataType the data type of the parameter
/// \param parameter_name name of the requested parameter
/// \param parameters list of parameters in which it will be searched
ParameterBase* findParameterByName(
    std::string const& parameter_name,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters);

/// Find an optional parameter of specific type for a given name.
///
/// \tparam ParameterDataType the data type of the parameter
/// \param parameter_name name of the requested parameter
/// \param parameters list of parameters in which it will be searched
/// \param num_components the number of components of the parameters or zero if
/// any number is acceptable
/// \param mesh an optional mesh pointer used for test whether the parameter is
/// defined on the given mesh. No test is performed if the pointer is a nullptr.
///
/// \see The documentation of the other findParameter() function.
template <typename ParameterDataType>
Parameter<ParameterDataType>* findParameterOptional(
    std::string const& parameter_name,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    int const num_components, MeshLib::Mesh const* const mesh = nullptr)
{
    // Find corresponding parameter by name.
    ParameterBase* parameter_ptr =
        findParameterByName(parameter_name, parameters);
    if (parameter_ptr == nullptr)
    {
        return nullptr;
    }

    // Check the type correctness of the found parameter.
    auto* const parameter =
        dynamic_cast<Parameter<ParameterDataType>*>(parameter_ptr);
    if (!parameter)
    {
        OGS_FATAL("The read parameter `{:s}' is of incompatible type.",
                  parameter_name);
    }

    if (num_components != 0 &&
        parameter->getNumberOfGlobalComponents() != num_components)
    {
        OGS_FATAL(
            "The read parameter `{:s}' has the wrong number of components "
            "({:d} instead of {:d}).",
            parameter_name, parameter->getNumberOfGlobalComponents(),
            num_components);
    }

    // Test the parameter's mesh only if there is a "test"-mesh provided.
    if (mesh != nullptr)
    {
        if (auto const error = isDefinedOnSameMesh(*parameter, *mesh))
        {
            OGS_FATAL(
                "The found parameter is not suitable for the use on the "
                "required mesh.\n{:s}",
                error->c_str());
        }
    }


    return parameter;
}

/// Find a parameter of specific type for a given name.
///
/// \tparam ParameterDataType the data type of the parameter
/// \param parameter_name name of the requested parameter
/// \param parameters list of parameters in which it will be searched
/// \param num_components the number of components of the parameters or zero if
/// any number is acceptable
/// \param mesh an optional mesh pointer used for test whether the parameter is
/// defined on the given mesh. No test is performed if the pointer is a nullptr.
///
/// \see The documentation of the other findParameter() function.
template <typename ParameterDataType>
Parameter<ParameterDataType>& findParameter(
    std::string const& parameter_name,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    int const num_components, MeshLib::Mesh const* const mesh = nullptr)
{
    auto* parameter = findParameterOptional<ParameterDataType>(
        parameter_name, parameters, num_components, mesh);

    if (!parameter)
    {
        OGS_FATAL(
            "Could not find parameter `{:s}' in the provided parameters list.",
            parameter_name);
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
    int const num_components, MeshLib::Mesh const* const mesh = nullptr)
{
    // Find parameter name in process config.
    //! \ogs_file_special
    auto const name = process_config.getConfigParameter<std::string>(tag);

    return findParameter<ParameterDataType>(name, parameters, num_components,
                                            mesh);
}

/// Find a parameter of specific type for a name given in the process
/// configuration under the optional tag.
/// If the tag is not present, a nullptr will be returned.
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
/// and return a pointer to that parameter. Additionally it checks for the
/// type of the found parameter.
template <typename ParameterDataType>
Parameter<ParameterDataType>* findOptionalTagParameter(
    BaseLib::ConfigTree const& process_config, std::string const& tag,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    int const num_components, MeshLib::Mesh const* const mesh = nullptr)
{
    // Find parameter name in process config.
    auto const name =
        //! \ogs_file_special
        process_config.getConfigParameterOptional<std::string>(tag);

    if (!name)
    {
        return nullptr;
    }
    return &findParameter<ParameterDataType>(*name, parameters, num_components,
                                             mesh);
}
}  // namespace ParameterLib
