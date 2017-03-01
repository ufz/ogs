/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "GroupBasedParameter.h"

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "MeshLib/Mesh.h"

namespace ProcessLib
{
std::unique_ptr<ParameterBase> createGroupBasedParameter(
    std::string const& name, BaseLib::ConfigTree const& config,
    MeshLib::Mesh const& mesh)
{
    //! \ogs_file_param{prj__parameters__parameter__type}
    config.checkConfigParameter("type", "Group");

    // get a property vector of group IDs
    //! \ogs_file_param{prj__parameters__parameter__Group__group_id_property}
    std::string const group_id_property_name = config.getConfigParameter<std::string>("group_id_property");
    DBUG("Using group_id_property %s", group_id_property_name.c_str());

    auto const& group_id_property =
        mesh.getProperties().getPropertyVector<int>(group_id_property_name);

    // parse mapping data
    typedef std::vector<double> Values;
    typedef std::pair<int, Values> Index_Values;
    std::vector<Index_Values> vec_index_values;
    //! \ogs_file_param{prj__parameters__parameter__Group__index_values}
    for (auto p : config.getConfigSubtreeList("index_values"))
    {
        //! \ogs_file_param{prj__parameters__parameter__Group__index_values__index}
        auto const index = p.getConfigParameter<int>("index");
        {
            //! \ogs_file_param{prj__parameters__parameter__Group__index_values__value}
            auto const value = p.getConfigParameterOptional<double>("value");

            if (value)
            {
                Values values(1, *value);
                vec_index_values.push_back(Index_Values(index, values));
                continue;
            }
        }

        // Value tag not available; continue with required values tag.
        //! \ogs_file_param{prj__parameters__parameter__Group__index_values__values}
        Values const values = p.getConfigParameter<Values>("values");

        if (values.empty())
            OGS_FATAL("No value available for constant parameter.");

        vec_index_values.push_back(Index_Values(index, values));
    }

    // check the input
    unsigned n_values = vec_index_values.front().second.size();
    for (auto p : vec_index_values)
    {
        auto itr = std::find(group_id_property->begin(), group_id_property->end(), p.first);
        if (itr == group_id_property->end())
            OGS_FATAL("Specified property index %d does not exist in the property vector %s",
                      p.first, group_id_property_name.c_str());

        if (p.second.size() != n_values)
            OGS_FATAL("The length of some values (%d) is different from the first one (%d). "
                      "The length should be same for all index_values.",
                      p.second.size(),  n_values);
    }

    // create a mapping table
    const int max_index =
        *std::max_element(group_id_property->begin(), group_id_property->end());
    std::vector<Values> vec_values(max_index + 1);
    for (auto p : vec_index_values)
        vec_values[p.first] = p.second;

    if (group_id_property->getMeshItemType() == MeshLib::MeshItemType::Node)
        return std::unique_ptr<ParameterBase>(
            new GroupBasedParameter<double, MeshLib::MeshItemType::Node>(
                name, *group_id_property, vec_values));
    else if (group_id_property->getMeshItemType() ==
             MeshLib::MeshItemType::Cell)
        return std::unique_ptr<ParameterBase>(
            new GroupBasedParameter<double, MeshLib::MeshItemType::Cell>(
                name, *group_id_property, vec_values));

    OGS_FATAL("Mesh item type of the specified property is not supported.");
}

}  // ProcessLib
