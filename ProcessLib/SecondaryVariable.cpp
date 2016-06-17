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
        boost::optional<BaseLib::ConfigTree> const& config,
        std::initializer_list<std::string> tag_names)
{
    if (!config) return;

    // read which variables are defined in the config
    for (auto const& tag_name : tag_names) {
        if (!_all_secondary_variables.insert(tag_name).second) {
            OGS_FATAL("Tag name <%s> has been specified twice as a secondary variable.");
        }

        //! \ogs_file_special
        if (auto var_name = config->getConfigParameterOptional<std::string>(tag_name))
        {
            // TODO check primary vars, too
            BaseLib::insertIfKeyValueUniqueElseError(
                        _map_tagname_to_varname, tag_name, *var_name,
                        "Secondary variable names must be unique.");
        }
    }
}

bool SecondaryVariableCollection::variableExists(std::string const& variable_name) const
{
    auto pred = [&variable_name](std::pair<std::string, std::string> const& p) {
        return p.second == variable_name;
    };

    // check if out_var is a  secondary variable
    auto const& var = std::find_if(
        _map_tagname_to_varname.cbegin(), _map_tagname_to_varname.cend(), pred);

    return var != _map_tagname_to_varname.cend();
}

void SecondaryVariableCollection::addSecondaryVariable(std::string const& tag_name,
                          const unsigned num_components,
                          SecondaryVariableFunctions&& fcts)
{
    auto it = _map_tagname_to_varname.find(tag_name);

    // get user-supplied var_name for the given tag_name
    if (it != _map_tagname_to_varname.end())
    {
        auto const& var_name = it->first;

        if (!_configured_secondary_variables
                 .emplace(std::make_pair(
                     var_name,
                     SecondaryVariable{
                         var_name, num_components, std::move(fcts)}))
                 .second)
        {
            OGS_FATAL("The secondary variable with name `%s' has already been "
                      "set up.",
                      var_name.c_str());
        }
    }
    else if (_all_secondary_variables.find(tag_name) ==
             _all_secondary_variables.end())
    {
        OGS_FATAL("The tag <%s> has not been registered to mark a secondary "
                  "variable.",
                  tag_name.c_str());
    }
}


} // namespace ProcessLib
