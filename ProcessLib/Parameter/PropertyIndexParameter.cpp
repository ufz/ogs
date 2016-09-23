/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "PropertyIndexParameter.h"

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "MeshLib/Mesh.h"

namespace ProcessLib
{
std::unique_ptr<ParameterBase> createPropertyIndexParameter(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& mesh)
{
    //! \ogs_file_param{parameter__type}
    config.checkConfigParameter("type", "PropertyIndex");

    // get a property vector
    //! \ogs_file_param{parameter_PropertyIndex__property_name}
    std::string const property_name = config.getConfigParameter<std::string>("property_name");
    DBUG("Using property_name %s", property_name.c_str());

    auto const& property_vector =
        mesh.getProperties().getPropertyVector<int>(property_name);
    if (!property_vector) {
        OGS_FATAL("The required property %s does not exist in the given mesh", property_name.c_str());
    }

    // parse mapping data
    typedef std::vector<double> Values;
    typedef std::pair<int, Values> Index_Values;
    std::vector<Index_Values> vec_index_values;
    //! \ogs_file_param{parameter_PropertyIndex__index_values}
    for (auto p : config.getConfigSubtreeList("index_values"))
    {
        //! \ogs_file_param{parameter_PropertyIndex__index_values__index}
        auto const index = p.getConfigParameter<int>("index");
        {
            //! \ogs_file_param{parameter_PropertyIndex__index_values__value}
            auto const value = p.getConfigParameterOptional<double>("value");

            if (value)
            {
                Values values(1, *value);
                vec_index_values.push_back(Index_Values(index, values));
                continue;
            }
        }

        // Value tag not available; continue with required values tag.
        //! \ogs_file_param{parameter_PropertyIndex__index_values__values}
        Values const values = p.getConfigParameter<Values>("values");

        if (values.empty())
            OGS_FATAL("No value available for constant parameter.");

        vec_index_values.push_back(Index_Values(index, values));
    }

    // check the input
    unsigned n_values = vec_index_values.front().second.size();
    for (auto p : vec_index_values)
    {
        auto itr = std::find(property_vector->begin(), property_vector->end(), p.first);
        if (itr == property_vector->end())
            OGS_FATAL("Specified property index %d does not exist in the property vector %s",
                      p.first, property_name.c_str());

        if (p.second.size() != n_values)
            OGS_FATAL("The length of some values (%d) is different from the first one (%d). "
                      "The length should be same for all index_values.",
                      p.second.size(),  n_values);
    }

    // create a mapping table
    const int max_index = *std::max_element(property_vector->begin(), property_vector->end());
    std::vector<Values> vec_values(max_index + 1);
    for (auto p : vec_index_values)
        vec_values[p.first] = p.second;

    if (property_vector->getMeshItemType() == MeshLib::MeshItemType::Node)
        return std::unique_ptr<ParameterBase>(
            new PropertyIndexParameter<double, MeshLib::MeshItemType::Node>(*property_vector, vec_values));
    else if (property_vector->getMeshItemType() == MeshLib::MeshItemType::Cell)
        return std::unique_ptr<ParameterBase>(
            new PropertyIndexParameter<double, MeshLib::MeshItemType::Cell>(*property_vector, vec_values));

    OGS_FATAL("Mesh item type of the specified property is not supported.");
}

}  // ProcessLib
