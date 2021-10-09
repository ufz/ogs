/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "GroupBasedParameter.h"

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "MeshLib/Mesh.h"

namespace ParameterLib
{
std::unique_ptr<ParameterBase> createGroupBasedParameter(
    std::string const& name, BaseLib::ConfigTree const& config,
    MeshLib::Mesh const& mesh)
{
    //! \ogs_file_param{prj__parameters__parameter__type}
    config.checkConfigParameter("type", "Group");

    // get a property vector of group IDs
    std::string const group_id_property_name =
        //! \ogs_file_param{prj__parameters__parameter__Group__group_id_property}
        config.getConfigParameter<std::string>("group_id_property");
    DBUG("Using group_id_property {:s}", group_id_property_name);

    auto const& group_id_property =
        mesh.getProperties().getPropertyVector<int>(group_id_property_name);

    // parse mapping data
    using Values = std::vector<double>;
    std::map<int, Values> vec_index_values;
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
                vec_index_values.emplace(index, values);
                continue;
            }
        }

        // Value tag not available; continue with required values tag.
        //! \ogs_file_param{prj__parameters__parameter__Group__index_values__values}
        Values const values = p.getConfigParameter<Values>("values");

        if (values.empty())
        {
            OGS_FATAL("No value available for constant parameter.");
        }

        vec_index_values.emplace(index, values);
    }

    // check the input
    for (auto p : vec_index_values)
    {
#ifndef USE_PETSC  // In case of partitioned meshes not all of the material ids
                   // might be available in the particular partition, therefore
                   // the check is omitted.
        if (std::find(group_id_property->begin(), group_id_property->end(),
                      p.first) == group_id_property->end())
        {
            OGS_FATAL(
                "Specified property index {:d} does not exist in the property "
                "vector {:s}",
                p.first, group_id_property_name);
        }
#endif

        auto const n_values = vec_index_values.begin()->second.size();
        if (p.second.size() != n_values)
        {
            OGS_FATAL(
                "The length of some values ({:d}) is different from the first "
                "one ({:d}). The length should be same for all index_values.",
                p.second.size(), n_values);
        }
    }

    if (group_id_property->getMeshItemType() == MeshLib::MeshItemType::Node)
    {
        return std::make_unique<
            GroupBasedParameter<double, MeshLib::MeshItemType::Node>>(
            name, mesh, *group_id_property, vec_index_values);
    }
    if (group_id_property->getMeshItemType() == MeshLib::MeshItemType::Cell)
    {
        return std::make_unique<
            GroupBasedParameter<double, MeshLib::MeshItemType::Cell>>(
            name, mesh, *group_id_property, vec_index_values);
    }

    OGS_FATAL("Mesh item type of the specified property is not supported.");
}

}  // namespace ParameterLib
