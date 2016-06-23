/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "SecondaryVariable.h"

namespace ProcessLib
{

SecondaryVariableCollection::SecondaryVariableCollection(
        std::initializer_list<std::string> internal_names)
{
    // read which variables are defined in the config
    for (auto const& internal_name : internal_names) {
        if (!_all_secondary_variables.insert(internal_name).second) {
            OGS_FATAL("Tag name <%s> has been specified twice as a secondary variable.");
        }
    }
}

void SecondaryVariableCollection::addNameMapping(std::string const& internal_name,
                                                 std::string const& external_name)
{
    _all_secondary_variables.insert(internal_name);

    // TODO check for missing secondary vars.
    // TODO check primary vars, too
    BaseLib::insertIfKeyUniqueElseError(
        _map_external_to_internal, external_name, internal_name,
        "Secondary variable names must be unique.");
}

bool SecondaryVariableCollection::variableExists(
    std::string const& external_name) const
{
    return _map_external_to_internal.find(external_name) !=
           _map_external_to_internal.cend();
}

void SecondaryVariableCollection::addSecondaryVariable(
    std::string const& internal_name,
    const unsigned num_components,
    SecondaryVariableFunctions&& fcts)
{
    if (!_configured_secondary_variables
             .emplace(std::make_pair(
                 internal_name,
                 SecondaryVariable{internal_name /* TODO change */,
                                   num_components, std::move(fcts)}))
             .second) {
        OGS_FATAL(
            "The secondary variable with internal name `%s' has already been "
            "set up.",
            internal_name.c_str());
    }
}

SecondaryVariable const& SecondaryVariableCollection::get(
    std::string const& external_name)
{
    auto const it = _map_external_to_internal.find(external_name);

    if (it == _map_external_to_internal.cend()) {
        OGS_FATAL(
            "A secondary variable with external name `%s' has not been set up.",
            external_name.c_str());
    }

    auto const& internal_name = it->second;
    auto const it2 = _configured_secondary_variables.find(internal_name);

    if (it2 == _configured_secondary_variables.end()) {

        OGS_FATAL(
            "A secondary variable with internal name `%s' has not been set up.",
            internal_name.c_str());
    }

    return it2->second;
}

} // namespace ProcessLib
