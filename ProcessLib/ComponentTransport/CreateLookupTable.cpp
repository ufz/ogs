/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateLookupTable.h"

#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <fstream>
#include <map>

#include "BaseLib/Algorithm.h"
#include "BaseLib/Error.h"
#include "BaseLib/FileTools.h"
#include "LookupTable.h"
#include "ProcessLib/ProcessVariable.h"

namespace
{
/**
 * returns indices of matching elements in a given vector.
 */
std::vector<std::size_t> getIndexVector(std::vector<double> const& data,
                                        double const value)
{
    std::vector<std::size_t> idx_vec;

    for (auto it = data.begin();
         (it = std::find(it, data.end(), value)) != data.end();
         ++it)
    {
        idx_vec.push_back(std::distance(data.begin(), it));
    }

    if (idx_vec.empty())
    {
        OGS_FATAL("No matching element {:d} is found in the vector.", value);
    }

    return idx_vec;
}

}  // namespace

namespace ProcessLib
{
namespace ComponentTransport
{
std::unique_ptr<LookupTable> createLookupTable(
    std::optional<std::string> const tabular_file,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>> const&
        process_variables)
{
    if (!tabular_file)
    {
        return nullptr;
    }

    auto const path_to_tabular_file =
        BaseLib::joinPaths(BaseLib::getProjectDirectory(), *tabular_file);

    if (!BaseLib::IsFileExisting(path_to_tabular_file))
    {
        OGS_FATAL(
            "Not found the tabular file with the specified file path: {:s}",
            path_to_tabular_file);
    }

    INFO("Found the tabular file: {:s}", path_to_tabular_file);

    std::ifstream in(path_to_tabular_file);
    if (!in)
    {
        OGS_FATAL("Couldn't open the tabular file: {:s}.",
                  path_to_tabular_file);
    }

    // read field names
    std::string line;
    std::getline(in, line);
    std::vector<std::string> field_names;
    boost::split(field_names, line, boost::is_any_of("\t "));

    // categorize entry fields
    std::vector<std::string> input_field_names;
    std::copy_if(field_names.begin(), field_names.end(),
                 std::back_inserter(input_field_names),
                 [](auto const& field_name) -> bool
                 { return field_name.find("_new") == std::string::npos; });

    // read table entries
    std::map<std::string, std::vector<double>> tabular_data;
    while (std::getline(in, line))
    {
        std::vector<std::string> field_data;
        boost::split(field_data, line, boost::is_any_of("\t "));

        assert(field_data.size() == field_names.size());
        for (std::size_t field_id = 0; field_id < field_data.size(); ++field_id)
        {
            tabular_data[field_names[field_id]].push_back(
                std::stod(field_data[field_id]));
        }
    }
    in.close();

    std::vector<Field> input_fields;
    input_fields.reserve(input_field_names.size());
    for (auto const& field_name : input_field_names)
    {
        auto pv = std::find_if(process_variables.begin(),
                               process_variables.end(),
                               [&field_name](auto const& v) -> bool
                               {
                                   return v[0].get().getName() == field_name ||
                                          v[0].get().getName() + "_prev" ==
                                              field_name;
                               });

        if (pv == process_variables.end())
        {
            OGS_FATAL(
                "Not found field name {:s} in the group of process variables "
                "defined in the project file.",
                field_name);
        }

        auto const process_id =
            static_cast<int>(std::distance(process_variables.cbegin(), pv));

        auto seed_points = tabular_data[field_name];
        BaseLib::makeVectorUnique(seed_points);

        std::vector<std::vector<std::size_t>> point_id_groups;
        for (auto const seed_point : seed_points)
        {
            auto const point_id_group =
                getIndexVector(tabular_data[field_name], seed_point);
            point_id_groups.push_back(point_id_group);
        }

        input_fields.emplace_back(point_id_groups, seed_points, field_name,
                                  process_id);
    }

    return std::make_unique<LookupTable>(std::move(input_fields),
                                         std::move(tabular_data));
}

}  // namespace ComponentTransport
}  // namespace ProcessLib
