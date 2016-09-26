/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MaterialIDParameter.h"

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/Properties.h"
#include "MeshLib/PropertyVector.h"

#include "ProcessLib/Parameter/MeshElementParameter.h"

namespace ProcessLib
{
std::unique_ptr<ParameterBase> createMaterialIDParameter(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& mesh)
{
    //! \ogs_file_param{parameter__type}
    config.checkConfigParameter("type", "MaterialID");

    // get a parameter name
    //! \ogs_file_param{parameter_MaterialID__name}
    std::string const parameter_name = config.getConfigParameter<std::string>("name");
    DBUG("Using parameter name %s", parameter_name.c_str());

    // get a property vector
    //! \ogs_file_param{parameter_MaterialID__property_name}
    std::string const index_property_name = config.getConfigParameter<std::string>("property_name");
    DBUG("Using property_name %s", index_property_name.c_str());

    auto const& index_property_vector =
        mesh.getProperties().getPropertyVector<int>(index_property_name);
    if (!index_property_vector) {
        OGS_FATAL("The required property %s does not exist in the given mesh",
                  index_property_name.c_str());
    }

    // parse mapping data
    typedef std::vector<double> Values;
    typedef std::pair<int, Values> Index_Values;
    std::vector<Index_Values> vec_index_values;
    //! \ogs_file_param{parameter_MaterialID__index_values}
    for (auto p : config.getConfigSubtreeList("index_values"))
    {
        //! \ogs_file_param{parameter_MaterialID__index_values__index}
        auto const index = p.getConfigParameter<int>("index");
        {
            //! \ogs_file_param{parameter_MaterialID__index_values__value}
            auto const value = p.getConfigParameterOptional<double>("value");

            if (value)
            {
                Values values(1, *value);
                vec_index_values.push_back(Index_Values(index, values));
                continue;
            }
        }

        // Value tag not available; continue with required values tag.
        //! \ogs_file_param{parameter_MaterialID__index_values__values}
        Values const values = p.getConfigParameter<Values>("values");

        if (values.empty())
            OGS_FATAL("No value available for constant parameter.");

        vec_index_values.push_back(Index_Values(index, values));
    }

    // check the input
    unsigned n_values = vec_index_values.front().second.size();
    for (auto p : vec_index_values)
    {
        auto itr = std::find(index_property_vector->begin(), index_property_vector->end(), p.first);
        if (itr == index_property_vector->end())
            OGS_FATAL("Specified property index %d does not exist in the property vector %s",
                      p.first, index_property_name.c_str());

        if (p.second.size() != n_values)
            OGS_FATAL("The length of some values (%d) is different from the first one (%d). "
                      "The length should be same for all index_values.",
                      p.second.size(),  n_values);
    }

    // create a mapping table
    const int max_index = *std::max_element(index_property_vector->begin(), index_property_vector->end());

    std::vector<std::size_t> prop_item2group_mapping(max_index+1);
    for (int i=0; i<max_index+1; i++)
        prop_item2group_mapping[i] = i;

    std::string prop_name = parameter_name;
    if (mesh.getProperties().hasPropertyVector(prop_name))
        prop_name += "_MaterialIDParameter";
    boost::optional<MeshLib::PropertyVector<double*> &> group_properties(
        const_cast<MeshLib::Mesh&>(mesh).getProperties().createNewPropertyVector<double*>(
            prop_name, prop_item2group_mapping.size(), prop_item2group_mapping,
            MeshLib::MeshItemType::Cell
        )
    );

    if (!group_properties)
        OGS_FATAL("Couldn't create a property vector in createMaterialIDParameter()");

    for (auto p : vec_index_values)
        (*group_properties).initPropertyValue(p.first, p.second);

    return std::unique_ptr<ParameterBase>(
        new MeshElementParameter<double*>(*group_properties));

}

}  // ProcessLib
